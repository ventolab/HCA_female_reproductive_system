"""
02_Run_Tacco.py

TACCO label transfer from postnatal single-cell reference to Visium HD data.

Outputs per donor in {VISIUM_HD_DIR}/{donor}/tacco_v2/:
  tacco_predictions.parquet        — barcode, coordinates, hard argmax celltype
  tacco_scores.parquet             — bins × cell types score matrix
  celltype_spatial_grid.pdf        — hard-label spatial grid
  celltype_score_spatial_grid.pdf  — TACCO score/probability spatial grid

Usage
-----
python 02_Run_Tacco.py
python 02_Run_Tacco.py --donor CW21
python 02_Run_Tacco.py --donor CW21 CW6
python 02_Run_Tacco.py --donor CW21 --overwrite
python 02_Run_Tacco.py --annotation-key fine_celltype
"""

import argparse
import math
import sys
from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import tacco as tc


VISIUM_HD_DIR   = Path("/lustre/scratch125/cellgen/vento/lm29/Tacco/data/VisiumHD/h5ad")
REFERENCE_PATH  = Path("/nfs/team292/projects/PanTissue/results/temp/01_integration/concat_postnatal_inner.h5ad")
ANNOTATIONS_CSV = Path("/nfs/team292/projects/PanTissue/results/freeze/annotations/concatenated_annotations_postnatal_v2.csv")

NCOLS           = 6
BG_COLOR        = "#DDDDDD"
BG_SIZE         = 1
FG_SIZE         = 3
DPI             = 200
SCORE_THRESHOLD = 0.4


# ── Reference loading ─────────────────────────────────────────────────────────

def load_reference(annotation_key: str) -> ad.AnnData:
    print(f"Loading reference: {REFERENCE_PATH} ...", flush=True)
    ref = sc.read(str(REFERENCE_PATH))
    print(f"  Reference shape: {ref.shape}", flush=True)

    print(f"Loading annotations: {ANNOTATIONS_CSV} ...", flush=True)
    ann = pd.read_csv(ANNOTATIONS_CSV, index_col="barcode")

    if annotation_key not in ann.columns:
        raise ValueError(
            f"annotation_key '{annotation_key}' not found in annotations CSV. "
            f"Available: {ann.columns.tolist()}"
        )

    if "cell_to_exclude" in ann.columns:
        ann = ann[~ann["cell_to_exclude"].astype(bool)]

    shared = ref.obs_names.intersection(ann.index)
    print(f"  Cells with annotations after exclusion: {len(shared):,} / {ref.n_obs:,}", flush=True)

    ref = ref[shared].copy()
    ref.obs[annotation_key] = ann.loc[shared, annotation_key].astype(str).values

    valid = ref.obs[annotation_key].str.lower() != "lowqc"
    ref   = ref[valid].copy()

    print(
        f"  Reference after filtering lowQC: {ref.n_obs:,} cells, "
        f"{ref.obs[annotation_key].nunique()} cell types",
        flush=True,
    )

    print(f"  Unique organs in reference: "
          f"{sorted(ref.obs['Organ'].astype(str).unique())}", flush=True)

    return ref


# ── Plotting ──────────────────────────────────────────────────────────────────

def _make_palette(celltypes: list) -> dict:
    cmap = plt.cm.get_cmap("tab20", 20)
    return {ct: cmap(i % 20) for i, ct in enumerate(sorted(celltypes))}


def plot_spatial_grid(xy: np.ndarray, cts: np.ndarray,
                      title: str, out_path: Path) -> None:
    unique_cts = sorted(set(cts))
    palette    = _make_palette(unique_cts)
    n          = len(unique_cts)
    nrows      = math.ceil(n / NCOLS)

    fig, axes = plt.subplots(nrows, NCOLS,
                             figsize=(6 * NCOLS, 5 * nrows),
                             squeeze=False)

    for idx, ct in enumerate(unique_cts):
        ax   = axes[idx // NCOLS][idx % NCOLS]
        mask = cts == ct
        ax.scatter(xy[~mask, 0], xy[~mask, 1],
                   s=BG_SIZE, c=BG_COLOR, rasterized=True)
        ax.scatter(xy[mask, 0],  xy[mask, 1],
                   s=FG_SIZE, c=[palette[ct]], rasterized=True)
        ax.set_title(ct, fontsize=10, pad=2)
        ax.axis("off")
        ax.invert_yaxis()

    for idx in range(n, nrows * NCOLS):
        axes[idx // NCOLS][idx % NCOLS].axis("off")

    fig.suptitle(title, fontsize=15, y=1.002)
    fig.tight_layout()
    fig.savefig(out_path, dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved hard-label plot → {out_path}", flush=True)


def plot_score_grid(xy: np.ndarray, scores: pd.DataFrame,
                    title: str, out_path: Path) -> None:
    celltypes = list(scores.columns)
    n         = len(celltypes)
    nrows     = math.ceil(n / NCOLS)

    fig, axes = plt.subplots(nrows, NCOLS,
                             figsize=(6 * NCOLS, 5 * nrows),
                             squeeze=False)

    for idx, ct in enumerate(celltypes):
        ax   = axes[idx // NCOLS][idx % NCOLS]
        vals = scores[ct].values
        sca  = ax.scatter(xy[:, 0], xy[:, 1],
                          c=vals, s=FG_SIZE, rasterized=True,
                          vmin=0, vmax=1)
        ax.set_title(ct, fontsize=10, pad=2)
        ax.axis("off")
        ax.invert_yaxis()
        fig.colorbar(sca, ax=ax, fraction=0.046, pad=0.04)

    for idx in range(n, nrows * NCOLS):
        axes[idx // NCOLS][idx % NCOLS].axis("off")

    fig.suptitle(title, fontsize=15, y=1.002)
    fig.tight_layout()
    fig.savefig(out_path, dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved score plot → {out_path}", flush=True)


# ── Per-donor annotation ──────────────────────────────────────────────────────

def annotate_one(query_path: Path, reference: ad.AnnData,
                 annotation_key: str, overwrite: bool) -> None:
    donor   = query_path.parent.name
    out_dir = query_path.parent / "tacco_v2"
    out_dir.mkdir(parents=True, exist_ok=True)

    preds_path      = out_dir / "tacco_predictions.parquet"
    scores_path     = out_dir / "tacco_scores.parquet"
    hard_plot_path  = out_dir / "celltype_spatial_grid.pdf"
    score_plot_path = out_dir / "celltype_score_spatial_grid.pdf"

    if preds_path.exists() and scores_path.exists() and score_plot_path.exists() and not overwrite:
        print(f"  Skipping {donor} — outputs already exist", flush=True)
        return

    print(f"\nLoading query: {query_path.name}  donor={donor}", flush=True)
    query = sc.read(str(query_path))
    print(f"  Query shape: {query.shape}", flush=True)

    if "counts" in query.layers:
        query.X = query.layers["counts"].copy()

    if "Organ" not in query.obs.columns:
        raise ValueError(f"Query {donor} does not contain query.obs['Organ'].")
    if "spatial" not in query.obsm:
        raise ValueError(f"Query {donor} does not contain query.obsm['spatial'].")

    organ      = str(query.obs["Organ"].iloc[0])
    print(f"  Query organ label: '{organ}'", flush=True)

    # case/whitespace-insensitive organ matching
    organ_norm  = organ.lower().replace("_", " ").strip()
    ref_organs  = reference.obs["Organ"].astype(str).str.lower().str.replace("_", " ").str.strip()
    reference_i = reference[ref_organs == organ_norm].copy()

    print(f"  Found {reference_i.n_obs:,} reference cells from Organ={organ}", flush=True)
    if reference_i.n_obs == 0:
        raise ValueError(
            f"No reference cells found for Organ='{organ}' (normalised: '{organ_norm}'). "
            f"Available organs in reference: {sorted(reference.obs['Organ'].astype(str).unique())}"
        )

    print("  Running tc.tl.annotate ...", flush=True)
    tc.tl.annotate(
        query,
        reference_i,
        annotation_key=annotation_key,
        result_key="predicted_celltype",
        bisections=0,
    )
    print("  Annotation complete.", flush=True)

    scores            = query.obsm["predicted_celltype"].copy()
    scores.index.name = "barcode"
    scores.columns    = scores.columns.astype(str)
    scores.index      = scores.index.astype(str)

    preds = pd.DataFrame({
        "barcode":            query.obs_names.astype(str),
        "pxl_col_in_fullres": query.obsm["spatial"][:, 0],
        "pxl_row_in_fullres": query.obsm["spatial"][:, 1],
        "celltype":           scores.idxmax(axis=1).values,
        "max_score":          scores.max(axis=1).values,
    })

    preds.to_parquet(preds_path, index=False)
    print(f"  Saved tacco_predictions.parquet ({len(preds):,} bins)", flush=True)

    scores.to_parquet(scores_path)
    print(f"  Saved tacco_scores.parquet "
          f"({scores.shape[0]:,} bins × {scores.shape[1]} cell types)", flush=True)

    xy         = query.obsm["spatial"]
    max_scores = scores.max(axis=1)
    high_conf  = max_scores > SCORE_THRESHOLD

    print(f"  High-confidence bins score>{SCORE_THRESHOLD}: "
          f"{high_conf.sum():,} / {len(high_conf):,}", flush=True)

    if high_conf.sum() > 0:
        plot_spatial_grid(
            xy[high_conf.values],
            preds.loc[high_conf.values, "celltype"].values,
            title=f"{donor} — TACCO {annotation_key} hard labels",
            out_path=hard_plot_path,
        )
    else:
        print("  WARNING: no high-confidence bins for hard-label plot", flush=True)

    plot_score_grid(
        xy, scores,
        title=f"{donor} — TACCO {annotation_key} scores",
        out_path=score_plot_path,
    )


# ── Main ──────────────────────────────────────────────────────────────────────

def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--annotation-key", default="fine_celltype",
                        help="Column from annotations CSV to use as cell type label.")
    parser.add_argument("--donor", default=None, nargs="+",
                        help="One or more donors (e.g. --donor CW21 CW6). If not set, runs all.")
    parser.add_argument("--overwrite", action="store_true",
                        help="Re-run even if outputs already exist.")
    args = parser.parse_args()

    reference = load_reference(args.annotation_key)

    if args.donor:
        query_paths = sorted(p for d in args.donor
                             for p in VISIUM_HD_DIR.glob(f"{d}/*_raw.h5ad"))
    else:
        query_paths = sorted(VISIUM_HD_DIR.glob("*/*_raw.h5ad"))

    if not query_paths:
        print(f"No *_raw.h5ad files found under {VISIUM_HD_DIR}", file=sys.stderr)
        sys.exit(1)

    print(f"\nFound {len(query_paths)} Visium HD sample(s):", flush=True)
    for p in query_paths:
        print(f"  {p.relative_to(VISIUM_HD_DIR)}", flush=True)

    for query_path in query_paths:
        try:
            annotate_one(
                query_path=query_path,
                reference=reference,
                annotation_key=args.annotation_key,
                overwrite=args.overwrite,
            )
        except Exception as e:
            import traceback
            print(f"\nERROR on {query_path.parent.name}: {type(e).__name__}: {e}",
                  file=sys.stderr)
            traceback.print_exc()

    print("\nDone.", flush=True)


if __name__ == "__main__":
    main()