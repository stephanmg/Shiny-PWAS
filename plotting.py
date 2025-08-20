import math

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots

from backend import get_single_gene_df, outcome_catalog

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
        if not values:
            continue  # nothing selected for this category -> don't constrain it

        kind = FILTER_MAP.get(ui_key)
        if not kind:
            continue  # unknown key

        is_kind = df["analysis_type"].eq(kind)

        if exact:
            # Exact match on Description
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

    if filters is None:
        return d
    else:
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


def plot_volcano(
    d: pd.DataFrame, use_log: bool, metric: str, limit: int, show_legend: bool
) -> go.Figure:
    fig = go.Figure()
    for gname, gdf in d.groupby("gene", sort=True):
        # Optional: nicer hover with extra columns if you have them
        hovertext = (
            "Gene: "
            + gdf["gene"].astype(str)
            + "<br>Type: "
            + gdf["analysis_type"].astype(str)
            + (
                ("<br>Desc: " + gdf["Description"].astype(str))
                if "Description" in gdf.columns
                else ""
            )
            + f"<br>{'−log10('+metric+')' if use_log else metric}: "
            + gdf["_y"].round(4).astype(str)
        )

        fig.add_trace(
            go.Scatter(
                x=gdf["_xj"],
                y=gdf["_y"],
                mode="markers",
                name=str(gname),
                hovertext=hovertext,
                hoverinfo="text",
                marker=dict(size=6, opacity=0.85),
            )
        )

    # Threshold line
    thresh_y = (-math.log10(0.05)) if use_log else 0.05
    fig.add_hline(y=thresh_y, line_dash="dot", line_width=1)

    # X axis ticks / labels from KIND_ORDER / KIND_LABEL
    fig.update_xaxes(
        tickmode="array",
        tickvals=list(range(len(KIND_ORDER))),
        ticktext=[KIND_LABEL[k] for k in KIND_ORDER],
    )

    # Titles, legend, layout
    fig.update_layout(
        title=f"{'−log10('+metric+')' if use_log else metric} by analysis type (top {limit}/gene/group)",
        xaxis_title="",
        yaxis_title=f"−log10({metric})" if use_log else metric,
        template="plotly_white",
        margin=dict(l=10, r=10, t=50, b=60),
        legend=dict(
            title="Gene",
            orientation="v",
            x=1.02,
            y=1.0,
            xanchor="left",
            yanchor="top",
            bgcolor="rgba(255,255,255,0.6)",
        ),
        showlegend=bool(show_legend),
    )
    return fig


def plot_bar(d: pd.DataFrame) -> go.Figure:
    # Count rows per (gene, analysis_type)
    counts = d.groupby(["gene", "analysis_type"]).size().reset_index(name="count")

    # Optional: ensure zero bars appear for missing (gene, category) combos
    # (uncomment if you want explicit zeros instead of missing bars)
    # genes = sorted(d["gene"].dropna().unique().tolist())
    # grid = pd.MultiIndex.from_product([genes, order],
    #                                   names=["gene", "analysis_type"]).to_frame(index=False)
    # counts = grid.merge(counts, on=["gene", "analysis_type"], how="left").fillna({"count": 0})

    counts["Category"] = (
        counts["analysis_type"].map(KIND_LABEL).fillna(counts["analysis_type"])
    )
    genes_order = sorted(counts["gene"].dropna().unique().tolist())
    cat_order = [KIND_LABEL[k] for k in KIND_ORDER]

    fig = px.bar(
        counts,
        x="gene",
        y="count",
        color="Category",
        barmode="group",
        category_orders={"gene": genes_order, "Category": cat_order},
        labels={"gene": "Gene", "count": "Number of outcomes"},
        title="Counts per gene by analysis type",
    )

    fig.update_layout(
        legend_title_text="Category",
        margin=dict(l=10, r=10, t=50, b=60),
    )
    return fig


def _to_float_series(s: pd.Series) -> pd.Series:
    return pd.to_numeric(s, errors="coerce")


def _safe_pthresh(pthresh):
    try:
        if pthresh in (None, "", "ALL"):
            return None
        return float(pthresh)
    except Exception:
        return None


def plot_heatmap(df, metric="p", use_log=True, pthresh=None) -> go.Figure:
    if df is None or df.empty:
        f = go.Figure()
        f.add_annotation(text="No data", x=0.5, y=0.5, showarrow=False)
        f.update_xaxes(visible=False)
        f.update_yaxes(visible=False)
        return f

    d = df.copy()
    if metric not in d.columns:
        metric = "p"

    m = pd.to_numeric(d[metric], errors="coerce")
    if pthresh not in (None, "", "ALL"):
        try:
            thr = float(pthresh)
            m = m.where(m <= thr)
        except Exception:
            pass

    m = m.clip(lower=1e-300, upper=1.0)
    d["_val"] = -np.log10(m) if use_log else m
    cbar_title = f"−log10({metric})" if use_log else metric

    fig = make_subplots(
        rows=2,
        cols=2,
        subplot_titles=[KIND_LABEL[k] for k in KIND_ORDER],
        horizontal_spacing=0.30,  # more space between columns
        vertical_spacing=0.35,  # more space between rows
    )

    pos = {
        KIND_ORDER[0]: (1, 1),
        KIND_ORDER[1]: (1, 2),
        KIND_ORDER[2]: (2, 1),
        KIND_ORDER[3]: (2, 2),
    }
    all_z = []

    for kind in KIND_ORDER:
        r, c = pos[kind]
        sub = d.loc[d["analysis_type"] == kind, ["gene", "Description", "_val"]].copy()
        if sub.empty:
            fig.add_annotation(
                text=f"No data for {KIND_LABEL[kind]}",
                x=0.5,
                y=0.5,
                xref=f"x{(r-1)*2+c}",
                yref=f"y{(r-1)*2+c}",
                showarrow=False,
            )
            continue

        mat = sub.pivot(index="Description", columns="gene", values="_val")

        z = mat.to_numpy(dtype=float)
        z[~np.isfinite(z)] = np.nan
        z_list = [[None if np.isnan(v) else float(v) for v in row] for row in z]
        all_z.extend([v for row in z_list for v in row if v is not None])

        xs = [str(x) for x in mat.columns.tolist()]
        ys = [str(y) for y in mat.index.tolist()]

        fig.add_trace(
            go.Heatmap(
                x=xs,
                y=ys,
                z=z_list,
                colorscale="Viridis",
                showscale=False,  # disable per-panel colorbar
                hovertemplate="Gene: %{x}<br>Phenotype: %{y}<br>"
                + cbar_title
                + ": %{z:.4f}"
                + f"<extra>{KIND_LABEL[kind]}</extra>",
            ),
            row=r,
            col=c,
        )
        fig.update_xaxes(tickangle=45, row=r, col=c)
        fig.update_yaxes(automargin=True, row=r, col=c)

    # Add ONE shared colorbar
    if all_z:
        fig.add_trace(
            go.Heatmap(
                z=[[min(all_z), max(all_z)]],  # dummy 2x1 matrix
                colorscale="Viridis",
                showscale=True,
                colorbar=dict(
                    title=cbar_title,
                    x=1.05,  # push it a bit outside the plots
                    y=0.5,
                    len=0.9,
                ),
                hoverinfo="skip",
            ),
            row=1,
            col=2,
        )

    fig.update_layout(template="plotly_white", margin=dict(l=40, r=100, t=80, b=40))
    return fig


def plot_bubble(d: pd.DataFrame, gene: str, category: str, metric: str) -> go.Figure:
    d = get_single_gene_df(d, gene, category)

    if d.empty:
        fig = go.Figure()
        fig.add_annotation(
            text=f"No rows for {gene} in {category}", x=0.5, y=0.5, showarrow=False
        )
        return fig

    if metric not in d.columns:
        metric = "p"  # fallback

    vals = d[metric].replace(0, 1e-300)
    d["neglog10"] = -np.log10(vals)

    fig = px.scatter(
        d,
        x="n_controls",
        y="Description",
        size="n_cases",
        color="neglog10",
        color_continuous_scale="Viridis_r",
        hover_data={"n_cases": True, "n_controls": True, metric: True},
        labels={
            "n_controls": "Number of controls",
            "Description": "Phenotype (Description)",
            "neglog10": f"−log10({metric})",
        },
        title=f"Bubble plot for {gene} (colored by −log10({metric}))",
    )

    fig.update_traces(marker=dict(line=dict(width=1, color="black")))
    fig.update_layout(
        yaxis=dict(tickfont=dict(size=8)),
        legend_title="n_cases",
        margin=dict(l=40, r=20, t=50, b=40),
    )
    return fig
