from unittest.mock import patch

import pandas as pd

import backend


class FakeResponse:
    """Class to mock requests"""

    def __init__(self, json_data, status=200):
        self._json = json_data
        self.status_code = status

    def raise_for_status(self):
        if self.status_code >= 400:
            raise RuntimeError(f"HTTP {self.status_code}")

    def json(self):
        return self._json


def fake_requests_get(url, timeout=30, params=None):
    """
    Minimal router for backend API calls used by backend.py.
    Adjust URLs as needed to match your backend.API base.
    """
    # Gene by Ensembl ID
    if "/gene/ensembl/ENSG000000001" in url:
        return FakeResponse({"ensembl_id": "ENSG000000001", "symbol": "GENE1"})

    # Gene by symbol/name (list response)
    if "/gene/name/PCSK9" in url:
        return FakeResponse([{"ensembl_id": "ENSGPCSK9", "symbol": "PCSK9"}])

    if "/gene/name/UNKNOWN" in url:
        return FakeResponse([])

    # Per-gene results
    if "/gene/ENSGPCSK9/results" in url:
        data = [
            # continuous
            {
                "analysis_type": "CONTINUOUS_VARIABLE",
                "outcome_id": "O1",
                "p": 0.01,
                "q": 0.02,
            },
            {
                "analysis_type": "CONTINUOUS_VARIABLE",
                "outcome_id": "O2",
                "p": 0.20,
                "q": 0.40,
            },
            # cv endpoints
            {
                "analysis_type": "CV_ENDPOINTS",
                "outcome_id": "O3",
                "p": 1e-10,
                "q": 5e-8,
            },
            # self reported
            {"analysis_type": "SELF_REPORTED", "outcome_id": "O4", "p": 0.03},
            # phecodes
            {"analysis_type": "PHECODES", "outcome_id": "O5", "q": 0.001},
        ]
        return FakeResponse({"results": data})

    # Outcome catalog
    if "/outcome" in url:
        cat = [
            {"id": "O1", "description": "HDL cholesterol"},
            {"id": "O2", "description": "LDL cholesterol"},
            {"id": "O3", "description": "Myocardial infarction"},
            {"id": "O4", "description": None, "name": "Self-reported X"},
            {"id": "O5", "outcome_string": "Phecode Y"},
        ]
        return FakeResponse(cat)

    return FakeResponse({"unknown": True}, status=404)


###############################################################################
# TESTS
###############################################################################


def test_parse_gene_list_basic():
    text = "METTL2A\n PCSK9,  APOB  \n\nENSG000000001"
    out = backend.parse_gene_list(text)
    assert out == ["METTL2A", "PCSK9", "APOB", "ENSG000000001"]


@patch("requests.get", side_effect=fake_requests_get)
def test_resolve_gene_symbol(mock_get):
    ensg, sym = backend.resolve_gene("PCSK9")
    assert ensg == "ENSGPCSK9"
    assert sym == "PCSK9"


@patch("requests.get", side_effect=fake_requests_get)
def test_resolve_gene_ensembl(mock_get):
    ensg, sym = backend.resolve_gene("ENSG000000001")
    assert ensg == "ENSG000000001"
    assert sym == "GENE1"


@patch("requests.get", side_effect=fake_requests_get)
def test_resolve_gene_unknown(mock_get):
    ensg, sym = backend.resolve_gene("UNKNOWN")
    assert ensg is None and sym is None


@patch("requests.get", side_effect=fake_requests_get)
def test_fetch_gene_results(mock_get):
    df = backend.fetch_gene_results("ENSGPCSK9", "BOTH")
    assert isinstance(df, pd.DataFrame)
    assert {"analysis_type", "outcome_id"}.issubset(df.columns)
    # sanity
    assert len(df) == 5


def _sample_results_df():
    return pd.DataFrame(
        [
            {
                "gene": "PCSK9",
                "analysis_type": "CONTINUOUS_VARIABLE",
                "outcome_id": "O1",
                "p": 0.01,
                "q": 0.02,
            },
            {
                "gene": "PCSK9",
                "analysis_type": "CONTINUOUS_VARIABLE",
                "outcome_id": "O2",
                "p": 0.20,
                "q": 0.40,
            },
            {
                "gene": "PCSK9",
                "analysis_type": "CV_ENDPOINTS",
                "outcome_id": "O3",
                "p": 1e-10,
                "q": 5e-8,
            },
            {
                "gene": "PCSK9",
                "analysis_type": "SELF_REPORTED",
                "outcome_id": "O4",
                "p": 0.03,
            },
            {
                "gene": "PCSK9",
                "analysis_type": "PHECODES",
                "outcome_id": "O5",
                "q": 0.001,
            },
            {
                "gene": "APOB",
                "analysis_type": "CONTINUOUS_VARIABLE",
                "outcome_id": "O1",
                "p": 0.05,
                "q": 0.06,
            },
            {
                "gene": "APOB",
                "analysis_type": "CV_ENDPOINTS",
                "outcome_id": "O3",
                "p": 0.9,
                "q": 0.9,
            },
        ]
    )


def _sample_catalog_df():
    return pd.DataFrame(
        [
            {"outcome_id": "O1", "description": "HDL cholesterol"},
            {"outcome_id": "O2", "description": "LDL cholesterol"},
            {"outcome_id": "O3", "description": "Myocardial infarction"},
            {"outcome_id": "O4", "name": "Self-reported X"},
            {"outcome_id": "O5", "outcome_string": "Phecode Y"},
        ]
    )


def test_enrich_labels_prefers_catalog_then_fallback(monkeypatch):
    df = _sample_results_df()

    # Stub outcome_catalog() to return our local frame (skip network + cache)
    monkeypatch.setattr(backend, "outcome_catalog", lambda: _sample_catalog_df())

    out = backend.enrich_labels(df)
    # ensure Description exists and is filled sensibly
    assert "Description" in out.columns
    # Known descriptions from catalog
    assert out.loc[out["outcome_id"] == "O1", "Description"].iat[0] == "HDL cholesterol"
    assert out.loc[out["outcome_id"] == "O2", "Description"].iat[0] == "LDL cholesterol"


def test_tidy_table_columns(monkeypatch):
    df = _sample_results_df()
    monkeypatch.setattr(backend, "outcome_catalog", lambda: _sample_catalog_df())

    tidy = backend.tidy_table(df)
    assert list(tidy.columns) == ["Gene", "Outcome ID", "Description", "p", "q"]
    # Same number of rows as input (tidy_table doesn’t limit)
    assert len(tidy) == len(df)
    # No NaN in Description (should be filled)
    assert tidy["Description"].fillna("").str.len().gt(0).all()


def test_top_n_per_gene_sorts_and_limits(monkeypatch):
    df = _sample_results_df()
    # choose continuous subset: 2 rows for PCSK9, 1 for APOB
    top1 = backend.top_n_per_gene(df, "CONTINUOUS_VARIABLE", n=1)
    # Should be at most 1 per gene
    counts = top1.groupby("gene").size().to_dict()
    assert counts.get("PCSK9", 0) <= 1
    assert counts.get("APOB", 0) <= 1

    # If we ask for 5, we still only get what's available per gene
    top5 = backend.top_n_per_gene(df, "CONTINUOUS_VARIABLE", n=5)
    counts5 = top5.groupby("gene").size().to_dict()
    assert counts5.get("PCSK9", 0) == 2
    assert counts5.get("APOB", 0) == 1


@patch("requests.get", side_effect=fake_requests_get)
def test_end_to_end_symbol_to_tidy_table(mock_get, monkeypatch):
    """
    Mini integration:
    - resolve PCSK9
    - fetch results
    - enrich/ tidy
    """
    ensg, sym = backend.resolve_gene("PCSK9")
    assert ensg == "ENSGPCSK9" and sym == "PCSK9"
    df = backend.fetch_gene_results(ensg, "BOTH")
    df["gene"] = sym
    # force catalog to our fake
    monkeypatch.setattr(backend, "outcome_catalog", lambda: _sample_catalog_df())

    tidy = backend.tidy_table(df)
    assert not tidy.empty
    assert "Description" in tidy.columns
    # sanity: at least the O3 MI row should be present
    assert (tidy["Outcome ID"] == "O3").any()
