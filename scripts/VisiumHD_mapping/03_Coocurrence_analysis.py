"""
03_Coocurrance_analysis.py

Soft co-occurrence analysis + colocalisation scatter + barplots.

Outputs in OUT_DIR/{organ_display}/{donor}/:
  coloc_{config_key}_{donor}.pdf/.png      — spatial scatter (hard labels, focal enlarged)
  cooccurrence_{config_key}_{donor}.pdf    — co-occurrence heatmap (all cell types)
  cooccurrence_barplot_{config_key}_{donor}.pdf  — focal vs neighbours bar chart
  cooccurrence_{config_key}_{donor}.csv    — long-format co-occurrence table

Usage
-----
python coloc_cooccurrence.py
python coloc_cooccurrence.py --donor CW21
python coloc_cooccurrence.py --organ Ovary
python coloc_cooccurrence.py --max-distance 100 --overwrite
"""

import argparse
import math
from pathlib import Path

import matplotlib
import matplotlib.lines
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
import seaborn as sns
import tacco as tc

matplotlib.rcParams["pdf.fonttype"] = 42

# ═══════════════════════════════════════════════════════════════════════════════
# ── USER CONFIG ───────────────────────────────────────────────────────────────
# ═══════════════════════════════════════════════════════════════════════════════

VISIUM_HD_DIR    = Path("/lustre/scratch125/cellgen/vento/lm29/Tacco/data/VisiumHD/h5ad")
OUT_DIR          = Path("/lustre/scratch125/cellgen/vento/lm29/Paper_plots/results/coloc_cooccurrence")

DPI              = 300
CAPTURE_AREA_UM  = 6500.0
SCORE_THRESHOLD  = 0.4
MIN_SPOTS        = 20
MAX_DISTANCE_DEFAULT = 100

SPOT_SIZE_BG     = 0.4
SPOT_SIZE_NB     = 1.5
SPOT_SIZE_FOCAL  = 10.0
SCALE_BAR_UM     = 500
SCALE_BAR_COLOR  = "black"

NCOLS            = 6
BG_COLOR         = "#DDDDDD"

VSMC_COLS        = ["Mesen_vSMCs_medArt", "Mesen_vSMCs_largeArt"]
COMPOSITION_KEY  = "predicted_celltype"

ORGAN_CONFIG = {
    "Ovary": {
        "organ_display": "Ovary",
        "donors": ["AU30", "AU14_AU18", "M23"],
        "focal":  "Mesen_AdvFibs_PI16low",
        "neighbours": [
            "vSMCs_combined",
            "Mesen_SMCs",
            "Mesen_OvarianFibs_Outcor",
            "Mesen_OvarianFibs_Perifol",
            "Mesen_OvarianFibs_InncorMedulla",
            "Meso_OSE",
        ],
        "colors": {
            "vSMCs_combined":                  "#2471A3",
            "Mesen_SMCs":                      "#F4D03F",
            "Mesen_OvarianFibs_Outcor":        "#8E44AD",
            "Mesen_OvarianFibs_Perifol":       "#F5CBA7",
            "Mesen_OvarianFibs_InncorMedulla": "#17A589",
            "Meso_OSE":                        "#90D5FF",
        },
        "focal_color": "#E74C3C",
    },
    "Fallopian_Tube_PI16low": {
        "organ_display": "Fallopian_Tube",
        "donors": ["CW21"],
        "focal":  "Mesen_AdvFibs_PI16low",
        "neighbours": [
            "vSMCs_combined",
            "Mesen_SMCs",
            "Epi_FallopianSec_EstrInd",
            "Mesen_FallopianFibs",
        ],
        "colors": {
            "vSMCs_combined":           "#2471A3",
            "Mesen_SMCs":               "#F4D03F",
            "Epi_FallopianSec_EstrInd": "#90D5FF",
            "Mesen_FallopianFibs":      "#17A589",
        },
        "focal_color": "#E74C3C",
    },
    "Fallopian_Tube_PI16hi": {
        "organ_display": "Fallopian_Tube",
        "donors": ["CW21"],
        "focal":  "Mesen_AdvFibs_PI16hi",
        "neighbours": [
            "vSMCs_combined",
            "Mesen_SMCs",
            "Epi_FallopianSec_EstrInd",
            "Mesen_FallopianFibs",
        ],
        "colors": {
            "vSMCs_combined":           "#2471A3",
            "Mesen_SMCs":               "#F4D03F",
            "Epi_FallopianSec_EstrInd": "#90D5FF",
            "Mesen_FallopianFibs":      "#17A589",
        },
        "focal_color": "#E74C3C",
    },
}

# ═══════════════════════════════════════════════════════════════════════════════
# ── END USER CONFIG ───────────────────────────────────────────────────────────
# ═══════════════════════════════════════════════════════════════════════════════


# ── CLI ───────────────────────────────────────────────────────────────────────

def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("--donor", default=None, nargs="+",
                    help="One or more donor IDs. Default: all.")
    ap.add_argument("--organ", default=None,
                    help="Run only this config key (e.g. Ovary). Default: all.")
    ap.add_argument("--max-distance", type=float, default=MAX_DISTANCE_DEFAULT,
                    help=f"Max co-occurrence distance in µm (default: {MAX_DISTANCE_DEFAULT}).")
    ap.add_argument("--overwrite", action="store_true")
    return ap.parse_args()


# ── Coordinate normalisation ──────────────────────────────────────────────────

def coords_to_um(xy: np.ndarray, donor: str) -> np.ndarray:
    max_range = float((xy.max(axis=0) - xy.min(axis=0)).max())
    if abs(max_range - CAPTURE_AREA_UM) / CAPTURE_AREA_UM < 0.5:
        return xy.copy()
    scale = CAPTURE_AREA_UM / max_range
    print(f"  [{donor}] scaled pixels→µm (scale={scale:.5f})", flush=True)
    return xy * scale


# ── Data loading + vSMC merge ─────────────────────────────────────────────────

def load_donor_scores(donor: str) -> tuple | None:
    """
    Load raw h5ad + tacco scores for one donor.
    Merges vSMC subtypes into vSMCs_combined BEFORE row-normalisation.
    Returns (xy_um, scores_df) or None if files are missing.
    """
    donor_dir   = VISIUM_HD_DIR / donor
    h5ad_path   = donor_dir / f"{donor}_raw.h5ad"
    scores_path = donor_dir / "tacco_v2" / "tacco_scores.parquet"

    for p in [h5ad_path, scores_path]:
        if not p.exists():
            print(f"  [{donor}] WARNING: {p.name} not found, skipping", flush=True)
            return None

    adata     = sc.read_h5ad(h5ad_path)
    scores_df = pd.read_parquet(scores_path)

    shared    = adata.obs_names.intersection(scores_df.index)
    adata     = adata[shared]
    scores_df = scores_df.loc[shared].copy()

    # ── merge vSMC subtypes → vSMCs_combined (before normalisation) ──────────
    vsmc_present = [c for c in VSMC_COLS if c in scores_df.columns]
    if vsmc_present:
        scores_df["vSMCs_combined"] = scores_df[vsmc_present].sum(axis=1)
        scores_df.drop(columns=vsmc_present, inplace=True)
        print(f"  [{donor}] merged {vsmc_present} → vSMCs_combined", flush=True)

    # row-normalise so scores sum to 1
    row_sums  = scores_df.sum(axis=1).replace(0, 1)
    scores_df = scores_df.div(row_sums, axis=0)

    xy_um = coords_to_um(adata.obsm["spatial"], donor)
    print(f"  [{donor}] {len(shared):,} bins, {scores_df.shape[1]} cell types", flush=True)
    return xy_um, scores_df


def build_anndata_for_cooccurrence(
    xy_um: np.ndarray,
    scores_df: pd.DataFrame,
    donor: str,
) -> ad.AnnData:
    """Wrap scores into an AnnData suitable for tc.tl.co_occurrence_matrix."""
    new_index          = scores_df.index.astype(str) + "_" + donor
    scores_df          = scores_df.copy()
    scores_df.index    = new_index

    result                      = ad.AnnData(X=np.ones((len(scores_df), 1), dtype=np.float32))
    result.obs_names            = new_index
    result.obsm[COMPOSITION_KEY] = scores_df
    result.obsm["spatial"]      = xy_um
    result.obs["donor"]         = pd.Categorical([donor] * len(result))
    return result


# ── Co-occurrence ─────────────────────────────────────────────────────────────

def compute_cooccurrence(adata: ad.AnnData, max_distance: float) -> pd.DataFrame:
    result  = tc.tl.co_occurrence_matrix(
        adata,
        annotation_key=COMPOSITION_KEY,
        position_key="spatial",
        max_distance=max_distance,
        sparse=True,
        reads=False,
        verbose=0,
    )
    occ_2d = result["occ"][:, :, 0]
    cts    = list(result["annotation"])
    return pd.DataFrame(occ_2d, index=cts, columns=cts)


def get_reliable_types(scores_df: pd.DataFrame) -> list:
    return [ct for ct in scores_df.columns if scores_df[ct].sum() >= MIN_SPOTS]


# ── Scale bar ─────────────────────────────────────────────────────────────────

def add_scale_bar(ax: plt.Axes, xy: np.ndarray,
                  bar_length_um: float = SCALE_BAR_UM,
                  color: str = SCALE_BAR_COLOR, fontsize: int = 7) -> None:
    x_min, x_max = float(xy[:, 0].min()), float(xy[:, 0].max())
    y_min, y_max = float(xy[:, 1].min()), float(xy[:, 1].max())
    pad_x = (x_max - x_min) * 0.03
    pad_y = (y_max - y_min) * 0.04
    bar_x_end, bar_x_start = x_max - pad_x, x_max - pad_x - bar_length_um
    bar_y = y_max - pad_y
    ax.plot([bar_x_start, bar_x_end], [bar_y, bar_y],
            color=color, linewidth=2.5, solid_capstyle="butt", zorder=10)
    ax.text((bar_x_start + bar_x_end) / 2, bar_y - pad_y * 0.5,
            f"{int(bar_length_um)} µm",
            ha="center", va="bottom", fontsize=fontsize, color=color, zorder=10)


# ── Colocalisation scatter ────────────────────────────────────────────────────

def plot_coloc_scatter(
    xy: np.ndarray,
    scores_df: pd.DataFrame,
    focal_key: str,
    focal_color: str,
    neighbour_keys: list,
    neighbour_colors: dict,
    donor: str,
    organ_display: str,
    config_key: str,
    out_path: Path,
) -> None:
    n           = len(xy)
    label       = np.full(n, "background", dtype=object)

    def get_score(key):
        if key in scores_df.columns:
            return scores_df[key].values.astype(np.float32)
        return np.zeros(n, dtype=np.float32)

    for nk in reversed(neighbour_keys):
        mask = get_score(nk) >= SCORE_THRESHOLD
        label[mask] = nk

    focal_mask = get_score(focal_key) >= SCORE_THRESHOLD
    label[focal_mask] = focal_key

    all_colors = {focal_key: focal_color, **neighbour_colors}
    color_arr  = np.full((n, 3), 0.88, dtype=np.float32)
    for key, hex_col in all_colors.items():
        m = label == key
        if m.any():
            rgb = np.array([int(hex_col.lstrip("#")[i:i+2], 16) / 255.0
                            for i in (0, 2, 4)], dtype=np.float32)
            color_arr[m] = rgb

    fig, (ax, ax_leg) = plt.subplots(
        1, 2, figsize=(9, 7),
        gridspec_kw={"width_ratios": [4, 1]},
    )

    bg = label == "background"
    if bg.any():
        ax.scatter(xy[bg, 0], xy[bg, 1], c=color_arr[bg],
                   s=SPOT_SIZE_BG, linewidths=0, rasterized=True, zorder=1)

    for nk in neighbour_keys:
        m = label == nk
        if m.any():
            ax.scatter(xy[m, 0], xy[m, 1], c=color_arr[m],
                       s=SPOT_SIZE_NB, linewidths=0, rasterized=True, zorder=2)

    if focal_mask.any():
        ax.scatter(xy[focal_mask, 0], xy[focal_mask, 1],
                   c=color_arr[focal_mask],
                   s=SPOT_SIZE_FOCAL, linewidths=0.2,
                   edgecolors="white", rasterized=True, zorder=3)

    add_scale_bar(ax, xy)
    ax.invert_yaxis()
    ax.axis("off")
    ax.set_title(
        f"{organ_display.replace('_', ' ')} — {donor}\nfocal: {focal_key}",
        fontsize=9,
    )

    ax_leg.axis("off")
    handles = [
        matplotlib.lines.Line2D([0], [0], marker="o", color="w",
                                 markerfacecolor=focal_color,
                                 markeredgecolor="white", markeredgewidth=0.2,
                                 markersize=11, label=focal_key),
    ]
    for nk in neighbour_keys:
        handles.append(matplotlib.lines.Line2D(
            [0], [0], marker="o", color="w",
            markerfacecolor=neighbour_colors[nk],
            markeredgewidth=0, markersize=7, label=nk,
        ))
    handles.append(matplotlib.lines.Line2D(
        [0], [0], marker="o", color="w",
        markerfacecolor="#E0E0E0", markeredgewidth=0, markersize=5,
        label=f"background (score<{SCORE_THRESHOLD})",
    ))
    ax_leg.legend(handles=handles, fontsize=7, framealpha=0.0,
                  loc="center left",
                  title=f"Score threshold: {SCORE_THRESHOLD}",
                  title_fontsize=7.5,
                  handletextpad=0.6, labelspacing=1.0, borderpad=0)

    fig.tight_layout()
    fig.savefig(out_path, bbox_inches="tight", dpi=DPI)
    fig.savefig(out_path.with_suffix(".png"), bbox_inches="tight", dpi=DPI)
    plt.close(fig)
    print(f"  Saved scatter: {out_path}", flush=True)


# ── Co-occurrence heatmap ─────────────────────────────────────────────────────

def plot_cooccurrence_heatmap(
    occ: pd.DataFrame,
    reliable_types: list,
    focal_key: str,
    neighbour_keys: list,
    title: str,
    out_path: Path,
) -> None:
    # restrict to focal + neighbours that exist in the occ matrix
    keep = [ct for ct in [focal_key] + neighbour_keys if ct in occ.index]
    if not keep:
        print("  WARNING: no focal/neighbour types in co-occurrence matrix, skipping heatmap.", flush=True)
        return

    sub      = occ.reindex(index=keep, columns=keep)
    log2_sub = np.log2(sub.clip(1e-10))

    unreliable = [ct for ct in keep if ct not in reliable_types]
    log2_sub.loc[unreliable] = np.nan
    log2_sub = log2_sub.clip(-3, 3)

    n_ct = len(keep)
    fig, ax = plt.subplots(figsize=(max(6, n_ct * 0.9), max(5, n_ct * 0.75)))
    cmap_hm = matplotlib.colormaps.get_cmap("RdBu_r").copy()
    cmap_hm.set_bad(color="#CCCCCC")

    vals  = log2_sub.values
    valid = vals[~np.isnan(vals)]
    vmax  = max(float(np.nanpercentile(np.abs(valid), 95)), 0.5) if len(valid) else 3.0

    sns.heatmap(log2_sub, ax=ax, cmap=cmap_hm, center=0,
                vmin=-vmax, vmax=vmax,
                linewidths=0.4, linecolor="white",
                cbar_kws={"label": "log₂(co-occurrence score)", "shrink": 0.6})
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right", fontsize=8)
    ax.set_yticklabels(ax.get_yticklabels(), rotation=0, fontsize=8)
    ax.set_title(title, fontsize=9)
    fig.tight_layout()
    fig.savefig(out_path, bbox_inches="tight", dpi=DPI)
    plt.close(fig)
    print(f"  Saved heatmap: {out_path}", flush=True)


# ── Co-occurrence barplot ─────────────────────────────────────────────────────

def plot_cooccurrence_barplot(
    occ: pd.DataFrame,
    reliable_types: list,
    focal_key: str,
    neighbour_keys: list,
    neighbour_colors: dict,
    focal_color: str,
    title: str,
    out_path: Path,
) -> None:
    """Bar chart of log2 co-occurrence scores: focal cell type vs each neighbour."""
    if focal_key not in occ.index:
        print(f"  WARNING: {focal_key} not in co-occurrence matrix, skipping barplot.", flush=True)
        return

    keep   = [ct for ct in neighbour_keys if ct in occ.columns]
    scores = np.log2(occ.loc[focal_key, keep].clip(1e-10))

    colors  = [neighbour_colors.get(ct, "#AAAAAA") for ct in keep]
    alphas  = [1.0 if ct in reliable_types else 0.35 for ct in keep]

    fig, ax = plt.subplots(figsize=(max(5, len(keep) * 0.9 + 1.5), 4))
    bars = ax.bar(range(len(keep)), scores, color=colors, edgecolor="white",
                  linewidth=0.5)
    for bar, alpha in zip(bars, alphas):
        bar.set_alpha(alpha)

    ax.axhline(0, color="black", linewidth=0.8, linestyle="--")
    ax.set_xticks(range(len(keep)))
    ax.set_xticklabels(keep, rotation=45, ha="right", fontsize=8)
    ax.set_ylabel("log₂(co-occurrence score)", fontsize=9)
    ax.set_title(title, fontsize=9)

    # legend for transparency
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor="gray", alpha=1.0, label=f"reliable (≥{MIN_SPOTS} eff. spots)"),
        Patch(facecolor="gray", alpha=0.35, label="low confidence"),
    ]
    ax.legend(handles=legend_elements, fontsize=7, framealpha=0.5)

    fig.tight_layout()
    fig.savefig(out_path, bbox_inches="tight", dpi=DPI)
    plt.close(fig)
    print(f"  Saved barplot: {out_path}", flush=True)


# ── CSV export ────────────────────────────────────────────────────────────────

def export_csv(occ: pd.DataFrame, focal_key: str,
               neighbour_keys: list, out_path: Path) -> None:
    keep = [ct for ct in [focal_key] + neighbour_keys if ct in occ.index]
    sub  = occ.reindex(index=keep, columns=keep)
    rows = [
        {
            "celltype1": ct1, "celltype2": ct2,
            "cooccurrence_score":      float(sub.loc[ct1, ct2]),
            "log2_cooccurrence_score": float(np.log2(max(sub.loc[ct1, ct2], 1e-10))),
        }
        for ct1 in keep for ct2 in keep
    ]
    (pd.DataFrame(rows)
       .sort_values("log2_cooccurrence_score", ascending=False)
       .to_csv(out_path, index=False))
    print(f"  Saved CSV: {out_path}", flush=True)


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    args          = parse_args()
    tag           = f"d{int(args.max_distance)}"
    configs_to_run = (
        {args.organ: ORGAN_CONFIG[args.organ]}
        if args.organ and args.organ in ORGAN_CONFIG
        else ORGAN_CONFIG
    )

    print(f"Max distance:  {args.max_distance} µm", flush=True)
    print(f"Score threshold: {SCORE_THRESHOLD}", flush=True)
    print(f"Output dir:    {OUT_DIR}", flush=True)

    for config_key, cfg in configs_to_run.items():
        organ_display    = cfg["organ_display"]
        donors           = args.donor if args.donor else cfg["donors"]
        focal_key        = cfg["focal"]
        focal_color      = cfg["focal_color"]
        neighbour_keys   = cfg["neighbours"]
        neighbour_colors = cfg["colors"]

        print(f"\n{'='*60}", flush=True)
        print(f"Config: {config_key}  focal: {focal_key}", flush=True)

        for donor in donors:
            print(f"\n  Donor: {donor}", flush=True)

            out_subdir = OUT_DIR / organ_display / donor
            out_subdir.mkdir(parents=True, exist_ok=True)

            scatter_path  = out_subdir / f"coloc_{config_key}_{donor}.pdf"
            heatmap_path  = out_subdir / f"cooccurrence_{config_key}_{tag}_{donor}.pdf"
            barplot_path  = out_subdir / f"cooccurrence_barplot_{config_key}_{tag}_{donor}.pdf"
            csv_path      = out_subdir / f"cooccurrence_{config_key}_{tag}_{donor}.csv"

            all_exist = all(p.exists() for p in [scatter_path, heatmap_path, barplot_path, csv_path])
            if all_exist and not args.overwrite:
                print(f"  [{donor}] all outputs exist, skipping. Use --overwrite to rerun.", flush=True)
                continue

            result = load_donor_scores(donor)
            if result is None:
                continue
            xy_um, scores_df = result

            reliable_types = get_reliable_types(scores_df)
            print(f"  [{donor}] reliable types (≥{MIN_SPOTS} eff. spots): {len(reliable_types)}", flush=True)

            # ── scatter plot ──────────────────────────────────────────────
            if not scatter_path.exists() or args.overwrite:
                plot_coloc_scatter(
                    xy_um, scores_df,
                    focal_key, focal_color,
                    neighbour_keys, neighbour_colors,
                    donor, organ_display, config_key,
                    scatter_path,
                )

            # ── co-occurrence ─────────────────────────────────────────────
            adata_occ = build_anndata_for_cooccurrence(xy_um, scores_df, donor)
            try:
                occ = compute_cooccurrence(adata_occ, args.max_distance)
            except Exception as e:
                print(f"  [{donor}] WARNING: co-occurrence failed: {e}", flush=True)
                continue

            export_csv(occ, focal_key, neighbour_keys, csv_path)

            title_base = (
                f"{organ_display.replace('_', ' ')} — {donor} — {config_key}\n"
                f"focal: {focal_key}  |  d ≤ {int(args.max_distance)} µm"
            )
            plot_cooccurrence_heatmap(
                occ, reliable_types,
                focal_key, neighbour_keys,
                title_base, heatmap_path,
            )
            plot_cooccurrence_barplot(
                occ, reliable_types,
                focal_key, neighbour_keys, neighbour_colors, focal_color,
                title_base, barplot_path,
            )

    print("\nDone.", flush=True)


if __name__ == "__main__":
    main()