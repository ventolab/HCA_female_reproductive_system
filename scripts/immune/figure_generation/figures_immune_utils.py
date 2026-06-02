import pandas as pd
import scanpy as sc
import numpy as np
import math
import matplotlib.pyplot as plt
from natsort import natsorted
from scipy.sparse import issparse
from scipy.stats import hypergeom
from statsmodels.stats.multitest import multipletests
import os
import seaborn as sns
from pathlib import Path


# Helper to filter markers present in adata and report status
def get_available_markers(marker_dict, name, adata):
    if not marker_dict:
        return {}
    
    avail = {
        ct: [g for g in genes if g in adata.var_names]
        for ct, genes in marker_dict.items()
    }
    # Remove empty categories
    avail = {ct: g for ct, g in avail.items() if g}
    
    all_genes = [g for genes in marker_dict.values() for g in genes]
    found_count = sum(len(v) for v in avail.values())
    missing = [g for g in all_genes if g not in adata.var_names]
    
    print(f"{name:20} found: {found_count} / {len(all_genes)}")
    if missing:
        print(f"  Not found: {missing}")
    return avail

# Bar plot
def barplot(
    which_var,
    adata,
    var              = "lineage",
    # --- figure geometry ---
    height           = 3,
    width            = 4,
    bar_height       = 0.8,         # thickness of each bar (0–1)
    # --- axes formatting ---
    xlabel           = "Proportion (%)",
    xtick_interval   = 10,          # spacing between x-axis ticks
    fontsize         = 8,
    # --- legend ---
    legend_cols      = 1,
    legend_x         = 1.02,        # bbox_to_anchor x
    legend_y         = 1.0,         # bbox_to_anchor y
    legend_fontsize  = 7,
    show_legend      = True,
    # --- output ---
    save_pdf         = False,
    save_png         = False,
    output_dir       = os.getcwd(),
    filename_prefix  = "",
    dpi              = 300,
):
    """
    Horizontal stacked bar chart showing percentage composition of `which_var`
    broken down by `var` (e.g. cell type proportions per tissue).

    Parameters
    ----------
    which_var       : str   Column in adata.obs that defines the stacked segments
                            (e.g. "celltype_HCA").
    adata           : AnnData
    var             : str   Column in adata.obs that defines the bars / rows
                            (e.g. "Tissue").
    height          : float Figure height in inches.
    width           : float Figure width in inches.
    bar_height      : float Bar thickness as a fraction of row spacing (0–1).
    xlabel          : str   X-axis label.
    xtick_interval  : int   Gap between x-axis tick marks.
    fontsize        : int   Font size for axis tick labels.
    legend_cols     : int   Number of columns in the legend.
    legend_x/y      : float Legend anchor position (in axes-fraction units).
    legend_fontsize : int   Font size for legend text.
    show_legend     : bool  Toggle legend visibility.
    save_pdf        : bool  Save figure as PDF.
    save_png        : bool  Save figure as PNG.
    output_dir      : str   Directory for saved files.
    filename_prefix : str   Optional prefix for the output filename.
    dpi             : int   Resolution for saved files.

    Returns
    -------
    fig, ax : matplotlib Figure and Axes objects for further customisation.
    """

    # ── 1. BUILD PERCENTAGE TABLE ──────────────────────────────────────────────
    # Cross-tabulate: rows = var categories, columns = which_var categories
    # normalize='index' makes each row sum to 100%
    plotdata = pd.crosstab(
        adata.obs[var],
        adata.obs[which_var],
        normalize="index"
    ) * 100

    # ── 2. RESPECT CATEGORY ORDER ─────────────────────────────────────────────
    # If var is a pandas Categorical, use its defined order (reversed so the
    # first category appears at the TOP of the horizontal bar chart).
    # Bug fix: the original forgot to assign the result of reorder_categories().
    if hasattr(adata.obs[var], "cat"):
        desired_order = adata.obs[var].cat.categories.tolist()[::-1]
        # Only reindex with categories that actually appear in the crosstab
        desired_order = [c for c in desired_order if c in plotdata.index]
        plotdata = plotdata.reindex(desired_order)

    # ── 3. RESOLVE COLOUR ORDER ───────────────────────────────────────────────
    # Pull colours from adata.uns
    # Build a dict {category_name: hex} first, then reindex to plotdata columns
    # so the colour order always matches the stacked segments regardless of how
    # the crosstab sorted its columns.
    uns_key = f"{which_var}_colors"

    if uns_key in adata.uns:
        categories   = adata.obs[which_var].cat.categories.tolist()
        uns_colors   = adata.uns[uns_key]                          # list of hex strings
        color_lookup = dict(zip(categories, uns_colors))           # {cell_type: hex}
        color_list   = [color_lookup.get(c, "#AAAAAA") for c in plotdata.columns]
    else:
        # Fallback: warn the user rather than failing silently
        print(f"WARNING: adata.uns['{uns_key}'] not found — plotting with matplotlib defaults. ")
        print(f"Run sc.pl.umap(adata, color='{which_var}') first to register colours, ")
        print(f"or set them manually with adata.uns['{uns_key}'].")
        color_list = None

    # ── 4. PLOT ───────────────────────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(width, height))

    plotdata.plot.barh(
        stacked    = True,
        ax         = ax,          # plot onto our axes (avoids figure duplication)
        edgecolor  = "none",
        color      = color_list,
        width      = bar_height,  # pandas barh uses 'width' for bar thickness
        legend     = False,       # we handle the legend ourselves below
        grid       = False
    )

    # ── 5. AXES FORMATTING ────────────────────────────────────────────────────
    ax.set_xticks(np.arange(0, 101, xtick_interval))
    ax.set_xlim(0, 100)
    ax.set_xlabel(xlabel, fontsize=fontsize)
    ax.set_ylabel(var, fontsize=fontsize)
    ax.tick_params(labelsize=fontsize)

    # Remove top and right spines for a cleaner publication look
    # ax.spines["top"].set_visible(False)
    # ax.spines["right"].set_visible(False)

    # ── 6. LEGEND ─────────────────────────────────────────────────────────────
    if show_legend:
        handles, labels = ax.get_legend_handles_labels()
        ax.legend(
            handles,
            labels,
            bbox_to_anchor  = (legend_x, legend_y),
            ncol            = legend_cols,
            fontsize        = legend_fontsize,
            frameon         = False,
            loc             = "upper left",
        )

    plt.tight_layout()

    # ── 7. SAVE ───────────────────────────────────────────────────────────────
    if save_pdf or save_png:
        safe_which = which_var.replace(" ", "_")
        safe_var   = var.replace(" ", "_")
        stem = os.path.join(
            output_dir,
            f"{filename_prefix}{safe_which}_per_{safe_var}"
        )
        if save_pdf:
            fig.savefig(stem + ".pdf", bbox_inches="tight", dpi=dpi)
        if save_png:
            fig.savefig(stem + ".png", bbox_inches="tight", dpi=dpi)

    return fig, ax