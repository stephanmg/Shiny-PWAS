# plotting.py
import numpy as np
import pandas as pd

KIND_ORDER = ["CONTINUOUS_VARIABLE","CV_ENDPOINTS","SELF_REPORTED","PHECODES"]
KIND_LABEL = {
    "CONTINUOUS_VARIABLE":"Continuous variables",
    "CV_ENDPOINTS":"Cardiovascular endpoints",
    "SELF_REPORTED":"Self reported",
    "PHECODES":"Phecodes",
}

def prepare_plot_df(df: pd.DataFrame, metric: str, limit: int) -> pd.DataFrame:
    """Return filtered/limited df with columns: analysis_type, gene, x, _xj, _y"""
    if df.empty or "analysis_type" not in df.columns: return pd.DataFrame()
    if "gene" not in df.columns: df = df.copy(); df["gene"] = "unknown"

    # valid metric rows
    d = df[df[metric].notna()].copy()
    # sort by q then p if present
    sort_cols = [c for c in ("q","p") if c in d.columns]
    d = d.sort_values(sort_cols or [metric], ascending=True)

    # take top N per (gene, analysis_type)
    d = d.groupby(["gene","analysis_type"], group_keys=False).head(max(1,int(limit)))
    d = d[d["analysis_type"].isin(KIND_ORDER)].copy()
    if d.empty: return d

    d["x"] = d["analysis_type"].map({k:i for i,k in enumerate(KIND_ORDER)}).astype(float)

    # safe values with epsilon for zeros when using -log10 in the app
    d["_vals"] = d[metric].astype(float).to_numpy()
    return d

def apply_scale(df: pd.DataFrame, use_log: bool, metric: str, eps: float=1e-300) -> pd.DataFrame:
    """Compute _y with optional -log10 and clip zeros."""
    if df.empty: return df
    vals = df["_vals"].copy()
    vals = np.clip(vals, eps, 1.0)
    if use_log:
        with np.errstate(divide="ignore", invalid="ignore"):
            df["_y"] = -np.log10(vals)
    else:
        df["_y"] = vals
    return df

def add_jitter(df: pd.DataFrame, seed: int=0, sd: float=0.045) -> pd.DataFrame:
    """ Adds jitter to plot """
    if df.empty: return df
    rng = np.random.default_rng(seed)
    df["_xj"] = df["x"] + rng.normal(0, sd, size=len(df))
    return df

