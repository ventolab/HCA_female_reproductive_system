"""
01_Prepare_VisiumHD.py

Loads Space Ranger outputs (square_08um binned) for each donor from a
metadata CSV, attaches spatial coordinates and sample metadata, applies
minimal QC filtering, and writes one <donor>_raw.h5ad per donor.

Output h5ad files are used as input by the TACCO annotation pipeline
(02_Run_Tacco.py).
"""

import argparse
from pathlib import Path
import pandas as pd
import scanpy as sc


def find_path_column(df):
    candidates = ["Location", "path", "spaceranger_path", "spaceranger_dir", "outs", "output_path"]
    for c in candidates:
        if c in df.columns:
            return c
    raise ValueError(
        f"Could not find path column. Expected one of: {candidates}\n"
        f"Found columns: {list(df.columns)}"
    )


def clean_columns(df):
    df = df.copy()
    df.columns = (
        df.columns.astype(str)
        .str.strip()
        .str.replace(" ", "_", regex=False)
        .str.replace("-", "_", regex=False)
    )

    if "Donor_id" in df.columns and "donor" not in df.columns:
        df = df.rename(columns={"Donor_id": "donor"})

    return df


def load_visium_hd_8um(row, path_col, outdir, overwrite=False):
    donor = str(row["donor"]).strip().replace(";", "_").replace(" ", "")

    sr_dir = Path(str(row[path_col]).strip().rstrip("/"))
    bin_dir = sr_dir / "binned_outputs" / "square_008um"
    matrix_h5 = bin_dir / "filtered_feature_bc_matrix.h5"
    spatial_dir = bin_dir / "spatial"

    outdir = Path(outdir) / donor
    out = outdir / f"{donor}_raw.h5ad"

    if out.exists() and not overwrite:
        print(f"Skipping existing file: {out}")
        return

    if not matrix_h5.exists():
        raise FileNotFoundError(f"Missing matrix file: {matrix_h5}")

    print(f"Loading donor={donor}")
    print(f"  {matrix_h5}")

    adata = sc.read_10x_h5(matrix_h5)
    adata.var_names_make_unique()

    # Add all CSV metadata columns to every bin
    for col, value in row.items():
        adata.obs[col] = value

    adata.obs["donor"] = donor
    adata.obs["bin_size"] = "square_008um"
    adata.obs["spaceranger_dir"] = str(sr_dir)

    # Add spatial coordinates
    parquet = spatial_dir / "tissue_positions.parquet"
    csv = spatial_dir / "tissue_positions.csv"

    if parquet.exists():
        pos = pd.read_parquet(parquet)
    elif csv.exists():
        pos = pd.read_csv(csv)
    else:
        pos = None
        print(f"  WARNING: no tissue_positions file found in {spatial_dir}")

    if pos is not None:
        pos = pos.set_index("barcode")
        adata.obs = adata.obs.join(pos, how="left")

        coord_cols = ["pxl_col_in_fullres", "pxl_row_in_fullres"]
        if all(c in adata.obs.columns for c in coord_cols):
            adata.obsm["spatial"] = adata.obs[coord_cols].to_numpy()
        else:
            print("  WARNING: pixel coordinate columns not found")

    image_files = (
        sorted(spatial_dir.glob("*.png"))
        + sorted(sr_dir.glob("*.tif"))
        + sorted(sr_dir.glob("*.tiff"))
    )

    adata.uns["spatial"] = {
        donor: {
            "metadata": {
                "donor": donor,
                "bin_size": "square_008um",
                "spaceranger_dir": str(sr_dir),
                "source_image_path": str(image_files[0]) if image_files else None,
            }
        }
    }

    # Minimal QC, keeping raw counts in .X
    sc.pp.filter_genes(adata, min_cells=3)
    sc.pp.filter_cells(adata, min_genes=100)
    sc.pp.calculate_qc_metrics(adata, inplace=True)

    outdir.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(out)

    print(f"  Saved: {out}")
    print(f"  Shape after QC: {adata.n_obs} bins x {adata.n_vars} genes")


def main():
    parser = argparse.ArgumentParser(
        description="Load Visium HD 8um Space Ranger outputs from CSV metadata and save h5ad files."
    )

    parser.add_argument(
        "--metadata",
        default="visiumhd_samples.csv",
        help="CSV file with sample metadata."
    )
    parser.add_argument(
        "--outdir",
        required=True,
        help="Output directory for h5ad files."
    )
    parser.add_argument(
        "--path-col",
        default=None,
        help="Column containing Space Ranger output paths. Defaults to Location."
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Overwrite existing h5ad files."
    )

    args = parser.parse_args()

    df = pd.read_csv(args.metadata)
    df = clean_columns(df)

    if "donor" not in df.columns:
        raise ValueError(f"Expected Donor_id/donor column. Found: {list(df.columns)}")

    path_col = args.path_col or find_path_column(df)
    outdir = Path(args.outdir)

    for _, row in df.iterrows():
        if row.isna().all():
            continue

        try:
            load_visium_hd_8um(
                row=row,
                path_col=path_col,
                outdir=outdir,
                overwrite=args.overwrite,
            )
        except Exception as e:
            print(f"ERROR processing donor={row.get('donor', 'UNKNOWN')}: {type(e).__name__}: {e}")


if __name__ == "__main__":
    main()