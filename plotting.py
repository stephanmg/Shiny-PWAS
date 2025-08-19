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
