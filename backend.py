from functools import lru_cache

import pandas as pd
import requests

API = "https://exphewas.statgen.org/v1/api"
CAT_KEEP = ["outcome_id", "description", "outcome_string", "name", "label", "phenotype"]


def parse_gene_list(text: str) -> list[str]:
    """Parse a list of genes, either separated by comma or newline
    Return: Gene symbol in upper case compatible for ExPheWas REST API
    """
    if not text:
        return []
    return [
        p.strip().upper() for p in text.replace(",", "\n").splitlines() if p.strip()
    ]


def resolve_gene(name: str) -> tuple[str | None, str | None]:
    """Return (ensg, symbol) or (None, None). Accepts symbol or ENSG."""
    try:
        if name.upper().startswith("ENSG"):
            r = requests.get(f"{API}/gene/ensembl/{name}", timeout=30)
            r.raise_for_status()
            j = r.json()
            return (j.get("ensembl_id") or name, j.get("symbol") or name)
        r = requests.get(f"{API}/gene/name/{name}", timeout=30)
        r.raise_for_status()
        j = r.json()
        if isinstance(j, list):
            if not j:
                return (None, None)
            j = j[0]
        return (j.get("ensembl_id"), j.get("symbol") or name)
    except Exception:
        return (None, None)


def fetch_gene_results(ensg: str, subset: str = "BOTH") -> pd.DataFrame:
    """Fetch results for a specified gene"""
    r = requests.get(
        f"{API}/gene/{ensg}/results", params={"analysis_subset": subset}, timeout=60
    )
    r.raise_for_status()
    data = r.json()
    rows = data["results"] if isinstance(data, dict) and "results" in data else data
    return pd.DataFrame(rows)


@lru_cache(maxsize=1)
def outcome_catalog() -> pd.DataFrame:
    """Retrieve outcome catalog for phenotyping categories"""
    r = requests.get(f"{API}/outcome", timeout=60)
    r.raise_for_status()
    cat = pd.DataFrame(r.json())
    if "id" in cat.columns:
        cat = cat.rename(columns={"id": "outcome_id"})
    keep = [c for c in CAT_KEEP if c in cat.columns]
    return cat[keep].copy()


def enrich_labels(df: pd.DataFrame) -> pd.DataFrame:
    """Left-join catalog & build a robust Description column."""
    if df.empty:
        return df.assign(Description=pd.Series(dtype="object"))
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
    return merged


def filter_table(df: pd.DataFrame, filters: list):
    """Filter table with a list of filter values"""
    if filters is None:
        return df
    if df is None or df.empty:
        return df

    mask = pd.Series(True, index=df.index)

    mask &= df["Description"].isin(filters)

    return df[mask].copy()


def tidy_table(
    df: pd.DataFrame,
    metric: str = "p",
    threshold: float = 0.05,
    filters: list = (),
) -> pd.DataFrame:
    """Format table"""
    if df.empty:
        return pd.DataFrame(columns=["Gene", "Outcome ID", "Description", "p", "q"])
    m = enrich_labels(df.copy())
    out = pd.DataFrame(
        {
            "Gene": m.get("gene", pd.Series([""] * len(m))),
            "Outcome ID": m["outcome_id"],
            "Description": m["Description"],
        }
    )
    if "p" in m.columns:
        out["p"] = m["p"]
    if "q" in m.columns:
        out["q"] = m["q"]

    if metric in out.columns:
        out = out[out[metric] < threshold]

    # if no filter provided, display all variables per default
    if filters is None:
        return out

    # if filter length is 0, meaning no variable selected to filter for, return identity
    if len(filters) == 0:
        return pd.DataFrame(columns=["Gene", "Outcome ID", "Description", "p", "q"])
    # else if filter not None, and variable provided by user, we filter for only these variables
    else:
        out = filter_table(out, filters)

    return out


def top_n_per_gene(df: pd.DataFrame, analysis_type: str, n: int) -> pd.DataFrame:
    """Top n samples / endpoints per gene"""
    sub = df[df["analysis_type"] == analysis_type].copy()
    if sub.empty:
        return sub
    sort_cols = [c for c in ("q", "p") if c in sub.columns]
    sub = sub.sort_values(sort_cols or ["p"], ascending=True)
    if "gene" not in sub.columns:
        sub["gene"] = "unknown"
    return sub.groupby("gene", group_keys=False).head(max(1, int(n)))


def get_single_gene_df(df: pd.DataFrame, gname: str, category: str) -> pd.DataFrame:
    """Return subset of df for a single gene chosen from the input list."""
    if df.empty or not gname:
        return pd.DataFrame()

    df = df[df["analysis_type"] == category]

    return df[df["gene"] == gname].copy()


###############################################################################
# CATALOG DATA
###############################################################################
KIND_CONT = "CONTINUOUS_VARIABLE"
KIND_CV = "CV_ENDPOINTS"
KIND_SELF = "SELF_REPORTED"
KIND_PHE = "PHECODES"
PATH_TO_CATALOG = "catalog/"

VALID_KINDS = {KIND_CONT, KIND_CV, KIND_SELF, KIND_PHE}


def _fetch_outcomes() -> pd.DataFrame:
    """Fetch the full outcomes catalog once."""
    r = requests.get(f"{API}/outcome", timeout=60)
    r.raise_for_status()
    return pd.DataFrame(r.json())


def _pick_label_column(df: pd.DataFrame) -> str:
    """
    Pick the best available human-readable label column.
    ExPheWAS often uses 'label'; some rows may have 'description'/'name'/'outcome_string'.
    """
    for col in ("label", "description", "name", "outcome_string"):
        if col in df.columns:
            return col
    # Fallback to 'id' if nothing else exists
    return "id"


def get_label_list(kind: str) -> list[str]:
    """
    Return a list of unique labels (Descriptions) for the requested analysis_type.
      kind ∈ {CONTINUOUS_VARIABLE, CV_ENDPOINTS, SELF_REPORTED, PHECODES}
    """
    if kind not in VALID_KINDS:
        raise ValueError(
            f"Unknown kind '{kind}'. Expected one of: {sorted(VALID_KINDS)}"
        )

    cat = _fetch_outcomes()
    if "analysis_type" not in cat.columns:
        return []

    sub = cat[cat["analysis_type"] == kind].copy()
    if sub.empty:
        return []

    label_col = _pick_label_column(sub)
    return (
        sub[label_col]
        .dropna()
        .astype(str)
        .str.strip()
        .replace({"": pd.NA})
        .dropna()
        .drop_duplicates()
        .sort_values(key=lambda s: s.str.lower())
        .tolist()
    )


###############################################################################
# HELPER FUNCTIONS
###############################################################################
def get_continuous_labels() -> list[str]:
    return get_label_list(KIND_CONT)


def get_cv_labels() -> list[str]:
    return get_label_list(KIND_CV)


def get_self_reported_labels() -> list[str]:
    return get_label_list(KIND_SELF)


def get_phecode_labels() -> list[str]:
    return get_label_list(KIND_PHE)


def get_all_label_lists() -> dict[str, list[str]]:
    return {
        KIND_CONT: get_continuous_labels(),
        KIND_CV: get_cv_labels(),
        KIND_SELF: get_self_reported_labels(),
        KIND_PHE: get_phecode_labels(),
    }
