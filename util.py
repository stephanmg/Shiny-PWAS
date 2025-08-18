import requests
from collections import defaultdict

API = "https://exphewas.statgen.org/v1/api"

def get_gene_associated_phenotypes(gene_name: str, analysis_subset: str = "BOTH"):
    """
    Returns a dict mapping analysis types to a list of (outcome_id, outcome_string).
    analysis_subset ∈ {"BOTH", "FEMALE_ONLY", "MALE_ONLY"}.
    """
    # 1) Resolve gene name -> Ensembl ID
    r = requests.get(f"{API}/gene/name/{gene_name}", timeout=30)
    r.raise_for_status()
    info = r.json()
    # The endpoint may return a list; pick the first match
    if isinstance(info, list):
        if not info:
            return {}
        info = info[0]
    ensg = info.get("ensembl_id")
    if not ensg:
        raise ValueError(f"Could not resolve Ensembl ID for '{gene_name}'")

    # 2) Fetch results for the gene
    r = requests.get(f"{API}/gene/{ensg}/results",
                     params={"analysis_subset": analysis_subset},
                     timeout=60)
    r.raise_for_status()
    data = r.json()
    # Some responses are wrapped under "results"
    rows = data["results"] if isinstance(data, dict) and "results" in data else data
    if not isinstance(rows, list):
        rows = []

    # 3) Group into types
    grouped = defaultdict(list)
    for row in rows:
        a_type = row.get("analysis_type")  # CONTINUOUS_VARIABLE, CV_ENDPOINTS, SELF_REPORTED, PHECODES
        oid = row.get("outcome_id")
        ostr = row.get("outcome_string") or row.get("description") or ""
        if a_type and oid:
            grouped[a_type].append((oid, ostr))

    # Sort each list by outcome_id for stability
    return {k: sorted(v, key=lambda t: str(t[0])) for k, v in grouped.items()}

if __name__ == "__main__":
    gene = "METTL2A"
    phenos = get_gene_associated_phenotypes(gene)
    print(phenos["CV_ENDPOINTS"][:5])
