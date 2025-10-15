import math

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from backend import outcome_catalog

KIND_ORDER = ["CONTINUOUS_VARIABLE", "CV_ENDPOINTS", "SELF_REPORTED", "PHECODES"]

KIND_LABEL = {
    "CONTINUOUS_VARIABLE": "Continuous variables",
    "CV_ENDPOINTS": "Cardiovascular endpoints",
    "SELF_REPORTED": "Self reported",
    "PHECODES": "Phecodes",
}

FILTER_MAP = {
    "filter_cont": "CONTINUOUS_VARIABLE",
    "filter_cv": "CV_ENDPOINTS",
    "filter_self": "SELF_REPORTED",
    "filter_phe": "PHECODES",
}


###############################################################################
# FILTER BY CATEGORY DESCRIPTION
###############################################################################
def filter_by_category_desc(
    df: pd.DataFrame,
    filters: dict,
    *,
    desc_col: str = "Description",
    exact: bool = True,
) -> pd.DataFrame:
    """
    Apply per-category filters on Description:
      - For rows where analysis_type == category, keep only those whose Description matches selected values
      - For other categories, leave rows unchanged
    `filters` should be like:
      {"filter_cont": [...], "filter_cv": [...], "filter_self": [...], "filter_phe": [...]}
    """
    if df is None or df.empty or not filters:
        return df
    if "analysis_type" not in df.columns or desc_col not in df.columns:
        return df

    mask = pd.Series(True, index=df.index)

    for ui_key, values in filters.items():
        if values is None:
            continue  # no filtering, display all entries, no checkbox ticked for filtering
        if len(values) == 0:
            pass

        kind = FILTER_MAP.get(ui_key)
        if not kind:
            continue  # unknown key

        is_kind = df["analysis_type"].eq(kind)

        if exact:
            # We match on empty string, so we get no matches, we invert this, because empty string means no variable selected for filtering
            is_match = df[desc_col].astype(str).isin(set(map(str, values)))
        else:
            # Case-insensitive substring match on Description
            s = df[desc_col].astype(str).str.lower()
            is_match = pd.Series(False, index=df.index)
            for v in values:
                is_match |= s.str.contains(str(v).lower(), na=False)

        # keep rows that are NOT this kind, OR (this kind AND matches)
        mask &= (~is_kind) | (is_kind & is_match)

    return df.loc[mask].copy()


###############################################################################
# HELPER FUNCTIONS FOR PLOTTING
###############################################################################
def prepare_plot_df(
    df: pd.DataFrame,
    metric: str,
    limit: int,
    catalog: pd.DataFrame,
    filters: dict = None,
) -> pd.DataFrame:
    """Return filtered/limited df with columns: analysis_type, gene, x, _xj, _y"""
    if df.empty or "analysis_type" not in df.columns:
        return pd.DataFrame()

    if "outcome_id" not in df.columns:
        if "outcome" in df.columns:
            df = df.rename(columns={"outcome": "outcome_id"})
        elif "id" in df.columns:
            df = df.rename(columns={"id": "outcome_id"})
    merged = df.merge(
        outcome_catalog(), on="outcome_id", how="left", suffixes=("", "_cat")
    )
    # pick first non-empty label
    for col in ("description", "outcome_string", "name", "label", "phenotype"):
        if col in merged.columns and merged[col].notna().any():
            desc = merged[col].fillna("")
            break
    else:
        desc = pd.Series([""] * len(merged))
    # fallback to results.outcome_string then outcome_id
    if "outcome_string" in df.columns:
        desc = desc.mask(desc.eq(""), df["outcome_string"])
    desc = desc.mask(desc.eq(""), merged["outcome_id"].astype(str))
    merged["Description"] = desc
    df = merged

    if "gene" not in df.columns:
        df = df.copy()
        df["gene"] = "unknown"

    # valid metric rows
    d = df[df[metric].notna()].copy()
    # sort by q then p if present
    sort_cols = [c for c in ("q", "p") if c in d.columns]
    d = d.sort_values(sort_cols or [metric], ascending=True)

    # take top N per (gene, analysis_type)
    d = d.groupby(["gene", "analysis_type"], group_keys=False).head(max(1, int(limit)))
    d = d[d["analysis_type"].isin(KIND_ORDER)].copy()

    if d.empty:
        return d

    d["x"] = (
        d["analysis_type"].map({k: i for i, k in enumerate(KIND_ORDER)}).astype(float)
    )

    # safe values with epsilon for zeros when using -log10 in the app
    d["_vals"] = d[metric].astype(float).to_numpy()

    return filter_by_category_desc(d, filters, desc_col="Description", exact=True)


def apply_scale(
    df: pd.DataFrame, use_log: bool, metric: str, eps: float = 1e-300
) -> pd.DataFrame:
    """Compute _y with optional -log10 and clip zeros."""
    if df.empty:
        return df
    vals = df["_vals"].copy()
    vals = np.clip(vals, eps, 1.0)
    if use_log:
        with np.errstate(divide="ignore", invalid="ignore"):
            df["_y"] = -np.log10(vals)
    else:
        df["_y"] = vals
    return df


def add_jitter(df: pd.DataFrame, seed: int = 0, sd: float = 0.045) -> pd.DataFrame:
    """Adds jitter to plot"""
    if df.empty:
        return df
    rng = np.random.default_rng(seed)
    df["_xj"] = df["x"] + rng.normal(0, sd, size=len(df))
    return df


###############################################################################
# PLOT TYPES FOR VISUALIZATION
###############################################################################


def heatmap_plot(
    df,
    metric="p",
    use_log=True,
    pthresh=None,
    nan_color="#FF00FF",
    single_plot=False,
    choice=None,
) -> matplotlib.figure.Figure:
    """
    Make 4 heatmaps (2x2), one per analysis_type.
    df must have columns: gene, analysis_type, Description, and the metric (p or q).
    """
    if df is None or df.empty:
        fig, ax = plt.subplots()
        ax.text(0.5, 0.5, "No data", ha="center", va="center")
        ax.set_axis_off()
        return fig

    d = df.copy()
    # filter threshold if given
    if pthresh is not None and metric in d.columns:
        d = d[d[metric] <= pthresh]

    if d.empty:
        fig, ax = plt.subplots()
        ax.text(0.5, 0.5, "No rows after filtering", ha="center", va="center")
        ax.set_axis_off()
        return fig

    # transform metric
    vals = d[metric].astype(float).to_numpy()
    vals = np.clip(vals, 1e-300, 1.0)  # avoid log10(0)
    d["_val"] = -np.log10(vals) if use_log else vals

    # draw grid
    if not single_plot:
        fig, axes = plt.subplots(2, 2, figsize=(12, 8))
        axes = axes.ravel()
    else:
        fig, axes = plt.subplots(1, 1, figsize=(12, 8))
        axes = [axes]

    vmin, vmax = d["_val"].min(), d["_val"].max()

    plots = KIND_ORDER
    if single_plot:
        if choice is not None:
            plots = [choice]

    for ax, kind in zip(axes, plots):
        sub = d[d["analysis_type"] == kind]
        print(sub.columns)
        if sub.empty:
            ax.text(
                0.5, 0.5, f"No data for\n{KIND_LABEL[kind]}", ha="center", va="center"
            )
            ax.set_axis_off()
            continue

        # Re-index to convert missing values to NaN values (otherwise matplotlib
        # renders missing values with default foreground color, which is white,
        # but we want to use the NaN color we set for the colormap via set_bad(...))
        x_axis_label = "outcome_id"  # was: Description (but too long)
        all_desc = sub[x_axis_label].unique()
        all_genes = sub["gene"].unique()
        mat = sub.pivot(index=x_axis_label, columns="gene", values="_val").reindex(
            index=all_desc, columns=all_genes
        )

        # set bad color to magenta (missing)
        cmap = plt.cm.viridis.copy()
        cmap.set_bad(color=nan_color)

        # plot finally
        im = ax.imshow(
            mat.values,
            aspect="auto",
            interpolation="nearest",
            vmin=vmin,
            vmax=vmax,
            cmap=cmap,
        )

        ax.set_title(KIND_LABEL[kind])
        ax.set_xticks(range(len(mat.columns)))
        ax.set_xticklabels(mat.columns, rotation=45, ha="right")
        ax.set_yticks(range(len(mat.index)))
        ax.set_yticklabels(mat.index, fontsize=4)

        cbar = fig.colorbar(im, ax=ax, shrink=0.7)
        cbar.set_label(f"-log10({metric})" if use_log else metric)

    return fig


def bubble_plot(
    d: pd.DataFrame,
    category: str,
    gene: str,
    metric: str,
    show_legend: bool,
) -> matplotlib.figure.Figure:
    """Creates a bubble plot for a single gene in a given category with provided metric"""
    if d.empty:
        fig, ax = plt.subplots()
        ax.text(0.5, 0.5, "No rows for single gene.", ha="center", va="center")
        ax.set_axis_off()
        return fig

    sizes = d["n_cases"] / d["n_cases"].max() * 800  # adjust scaling factor

    # CONTINUOUS_VARIABLE do not have size
    no_sizes = False
    if category == "CONTINUOUS_VARIABLE":
        no_sizes = True

    if metric not in d.columns:
        metric = "p"  # fallback
    vals = d[metric].replace(0, 1e-300)  # avoid log(0)

    # Use -log10 for better spread
    colors = -np.log10(vals)

    x = "n_controls"
    if category == "CONTINUOUS_VARIABLE":
        x = "n"

    print(d.columns)

    fig, ax = plt.subplots(figsize=(10, 5))
    sc = ax.scatter(
        d[x],
        d["Description"],
        s=sizes if not no_sizes else None,
        c=colors,
        cmap="viridis_r",
        alpha=0.8,
        edgecolor="k",
    )

    print("there")

    ax.set_xlabel("Number of controls")
    ax.set_ylabel("Phenotype (Description)")
    ax.set_title(f"Bubble plot for {gene} (colored by −log10({metric}))")
    ax.grid(True, linestyle="--", alpha=0.5)
    ax.tick_params(axis="y", labelsize=8)

    # Colorbar for metric
    cbar = fig.colorbar(sc, ax=ax, shrink=0.8)
    cbar.set_label(f"−log10({metric})")

    print("category:")
    print(category)
    # Legend for bubble sizes
    if show_legend and not no_sizes:
        for size in [d["n_cases"].min(), d["n_cases"].median(), d["n_cases"].max()]:
            ax.scatter(
                [],
                [],
                s=size / d["n_cases"].max() * 800,
                c="gray",
                alpha=0.5,
                edgecolor="k",
                label=f"{int(size)} cases",
            )
        ax.legend(scatterpoints=1, frameon=True, labelspacing=1, title="n_cases")

    fig.tight_layout()
    return fig


def bar_plot(d: pd.DataFrame, show_legend: bool) -> matplotlib.figure.Figure:
    """Creates 4 bar plots for all available phenotyping categories"""
    counts = d.groupby(["gene", "analysis_type"]).size().reset_index(name="count")

    genes = sorted(d["gene"].unique())
    x = np.arange(len(genes))
    width = 0.2

    fig, ax = plt.subplots(figsize=(8, 4.5))

    for i, kind in enumerate(KIND_ORDER):
        subset = counts[counts["analysis_type"] == kind]
        # Align with all genes (fill missing with 0)
        y = [subset.loc[subset["gene"] == g, "count"].sum() for g in genes]
        ax.bar(x + i * width, y, width, label=KIND_LABEL[kind])

    ax.set_xticks(x + width * (len(KIND_ORDER) - 1) / 2)
    ax.set_xticklabels(genes, rotation=45, ha="right")
    ax.set_ylabel("Number of outcomes")
    ax.set_title("Counts per gene by analysis type")
    if show_legend:
        ax.legend(title="Category")
    fig.tight_layout()
    return fig


def volcano_plot(
    d: pd.DataFrame,
    limit: int,
    use_log: bool = False,
    show_legend: bool = True,
    metric: str = "p",
) -> matplotlib.figure.Figure:
    """Create a grouped volcano plot for all available phenotyping categories"""
    fig, ax = plt.subplots(figsize=(7.2, 3.8))
    for gname, gdf in d.groupby("gene", sort=True):
        ax.scatter(gdf["_xj"], gdf["_y"], s=20, alpha=0.8, label=gname)

    thresh_y = -math.log10(0.05) if use_log else 0.05
    ax.axhline(thresh_y, linestyle=":", linewidth=1)
    ax.set_xticks(range(len(KIND_ORDER)))
    ax.set_xticklabels([KIND_LABEL[k] for k in KIND_ORDER])
    ax.set_ylabel(f"−log10({metric})" if use_log else metric)
    ax.set_title(
        f"{'−log10('+metric+')' if use_log else metric} by analysis type (top {limit}/gene/group)"
    )
    if show_legend:
        ax.legend(title="Gene", loc="upper right", frameon=True, fontsize="small")
    fig.tight_layout()
    return fig
