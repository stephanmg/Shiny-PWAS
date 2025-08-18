from shiny import App, ui, render, reactive
import requests
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from time import strftime
import math

API = "https://exphewas.statgen.org/v1/api"

app_ui = ui.page_sidebar(
    ui.sidebar(
        ui.input_text_area(
            "genes", "Gene names (comma or newline separated)",
            value="METTL2A\nPCSK9", width="100%", height="120px"
        ),
        ui.input_radio_buttons(
            "subset", "Analysis subset",
            choices={"BOTH": "Both", "MALE_ONLY": "Male only", "FEMALE_ONLY": "Female only"},
            selected="BOTH",
        ),
        ui.input_numeric("limit", "Max rows per gene per group", value=10, min=1, step=1),
        ui.input_action_button("btn_phenos", "Load phenotypes"),
        ui.input_select(
            "metric", "Value to plot",
            choices={"p": "p-value", "q": "q-value (FDR)"},
            selected="p",
        ),
        ui.input_checkbox("neglog10", "Use −log10 scale", value=True),
        ui.input_checkbox("show_legend", "Show legend", value=True),
        ui.input_action_button("btn_plot", "Load plot"),
        width=340,
    ),
    ui.tags.style(""".content-wrapper {max-width: 1100px; margin: 0 auto;}"""),
    ui.div(
        {"class": "content-wrapper"},
        ui.card(
            ui.card_header("Plot of selected metric by analysis type (color-coded by gene)"),
            ui.output_plot("plot_out", height="380px"),
        ),
        ui.card(
            ui.card_header("Phenotypes (top N per gene for each endpoint)"),
            ui.navset_pill(
                ui.nav_panel("Continuous variables", ui.output_table("tbl_continuous")),
                ui.nav_panel("Cardiovascular endpoints", ui.output_table("tbl_cv")),
                ui.nav_panel("Self reported", ui.output_table("tbl_self")),
                ui.nav_panel("Phecodes", ui.output_table("tbl_phecode")),
            ),
        ),
        ui.card(
            ui.card_header("Download combined table (ALL rows)"),
            ui.input_action_button("btn_download_csv", "Download CSV"),
            ui.output_text_verbatim("download_msg"),
        ),
        ui.card(
            ui.card_header("Log"),
            ui.output_text_verbatim("log_text"),
        ),
    ),
    title="ExPheWAS: multi-gene phenotypes view",
)

def log_append(logs: reactive.Value, msg: str):
    logs.set((logs.get() or "") + f"[{strftime('%H:%M:%S')}] {msg}\n")

def parse_gene_list(text: str) -> list[str]:
    if not text:
        return []
    parts = [p.strip() for p in text.replace(",", "\n").splitlines()]
    return [p for p in parts if p]

def resolve_gene(name: str):
    try:
        if name.upper().startswith("ENSG"):
            r = requests.get(f"{API}/gene/ensembl/{name}", timeout=30)
            r.raise_for_status()
            info = r.json()
            ensg = info.get("ensembl_id") or name
            symbol = info.get("symbol") or name
            return ensg, symbol
        else:
            r = requests.get(f"{API}/gene/name/{name}", timeout=30)
            r.raise_for_status()
            info = r.json()
            if isinstance(info, list):
                if not info:
                    return None, None
                info = info[0]
            ensg = info.get("ensembl_id")
            symbol = info.get("symbol") or name
            return ensg, symbol
    except Exception:
        return None, None

def fetch_gene_results(ensg: str, subset: str) -> pd.DataFrame:
    r = requests.get(f"{API}/gene/{ensg}/results", params={"analysis_subset": subset}, timeout=60)
    r.raise_for_status()
    data = r.json()
    rows = data["results"] if isinstance(data, dict) and "results" in data else data
    return pd.DataFrame(rows)

def fetch_outcome_catalog() -> pd.DataFrame:
    r = requests.get(f"{API}/outcome", timeout=60)
    r.raise_for_status()
    cat = pd.DataFrame(r.json())
    if "id" in cat.columns:
        cat = cat.rename(columns={"id": "outcome_id"})
    keep = ["outcome_id", "description", "outcome_string", "name", "label", "phenotype"]
    keep = [c for c in keep if c in cat.columns]
    return cat[keep]

OUTCOME_CATALOG = fetch_outcome_catalog()

def choose_label(df: pd.DataFrame, col_candidates=("description","outcome_string","name","label","phenotype")):
    for c in col_candidates:
        if c in df.columns:
            s = df[c]
            if s.notna().any():
                return s.fillna("")
    return pd.Series([""] * len(df))

def tidy_for_table(df: pd.DataFrame) -> pd.DataFrame:
    if df.empty:
        return pd.DataFrame(columns=["Gene", "Outcome ID", "Description", "p", "q"])
    if "outcome_id" not in df.columns:
        if "outcome" in df.columns:
            df = df.rename(columns={"outcome": "outcome_id"})
        elif "id" in df.columns:
            df = df.rename(columns={"id": "outcome_id"})
    merged = df.merge(OUTCOME_CATALOG, on="outcome_id", how="left", suffixes=("", "_cat"))
    desc = choose_label(merged, ("description", "outcome_string", "name", "label", "phenotype"))
    if "outcome_string" in df.columns:
        desc = desc.mask(desc.eq(""), df["outcome_string"])
    desc = desc.mask(desc.eq(""), merged["outcome_id"].astype(str))
    out = pd.DataFrame({
        "Gene": merged.get("gene", pd.Series([""] * len(merged))),
        "Outcome ID": merged["outcome_id"],
        "Description": desc,
    })
    if "p" in merged.columns:
        out["p"] = merged["p"]
    if "q" in merged.columns:
        out["q"] = merged["q"]
    return out

def server(input, output, session):
    logs = reactive.Value("")
    df_results = reactive.Value(pd.DataFrame())
    dl_msg = reactive.Value("Click 'Download CSV' to save the combined table with ALL rows.")

    @reactive.Effect
    @reactive.event(input.btn_phenos)
    def load_phenos():
        text = input.genes() or ""
        subset = input.subset()
        genes = parse_gene_list(text)
        if not genes:
            log_append(logs, "Enter at least one gene (symbol or ENSG).")
            df_results.set(pd.DataFrame())
            return
        frames = []
        for g in genes:
            log_append(logs, f"Resolving '{g}' → Ensembl ID…")
            ensg, symbol = resolve_gene(g)
            if not ensg:
                log_append(logs, f"  ! Could not resolve '{g}' — skipping.")
                continue
            try:
                log_append(logs, f"Fetching associations for {symbol} ({ensg}) …")
                df = fetch_gene_results(ensg, subset)
                if not df.empty:
                    df = df.copy()
                    df["gene"] = symbol
                    frames.append(df)
                    log_append(logs, f"  {len(df)} rows.")
                else:
                    log_append(logs, f"  0 rows.")
            except Exception as e:
                log_append(logs, f"  Error fetching for {symbol}: {e}")
        if frames:
            all_df = pd.concat(frames, ignore_index=True)
            df_results.set(all_df)
            log_append(logs, f"Total rows combined: {len(all_df)}")
        else:
            df_results.set(pd.DataFrame())
            log_append(logs, "No data loaded.")

    @reactive.Effect
    @reactive.event(input.btn_download_csv)
    def _download_all():
        df = df_results.get()
        if df is None or df.empty:
            dl_msg.set("No data to download.")
            return
        full = tidy_for_table(df.copy())
        path = "/mnt/data/multigene_results_all.csv"
        full.to_csv(path, index=False)
        dl_msg.set(f"Saved CSV to {path}")

    @output
    @render.text
    def download_msg():
        return dl_msg.get()

    @reactive.Effect
    @reactive.event(input.btn_plot)
    def _log_plot_click():
        log_append(logs, "Plot requested.")

    @output
    @render.plot
    def plot_out():
        df = df_results.get()
        if df is None or df.empty:
            fig, ax = plt.subplots()
            ax.text(0.5, 0.5, "No data yet. Click 'Load phenotypes'.", ha="center", va="center")
            ax.set_axis_off()
            return fig

        metric = input.metric()
        use_log = bool(input.neglog10())
        limit = max(1, int(input.limit()))

        if metric not in df.columns:
            fig, ax = plt.subplots()
            ax.text(0.5, 0.5, f"Column '{metric}' not found in results.", ha="center", va="center")
            ax.set_axis_off()
            return fig

        if "gene" not in df.columns:
            df = df.copy()
            df["gene"] = "unknown"

        plot_df = df.copy()
        plot_df = plot_df[pd.notnull(plot_df[metric])].copy()

        sort_cols = [c for c in ["q", "p"] if c in plot_df.columns]
        if sort_cols:
            plot_df = plot_df.sort_values(sort_cols, ascending=[True] * len(sort_cols))
        else:
            plot_df = plot_df.sort_values(metric, ascending=True)

        if "analysis_type" not in plot_df.columns:
            fig, ax = plt.subplots()
            ax.text(0.5, 0.5, "No 'analysis_type' column present.", ha="center", va="center")
            ax.set_axis_off()
            return fig

        # Limit per (gene, analysis_type)
        plot_df = plot_df.groupby(["gene", "analysis_type"], group_keys=False).head(limit)

        # Reduce to recognized types and assign x positions
        kind_order = ["CONTINUOUS_VARIABLE", "CV_ENDPOINTS", "SELF_REPORTED", "PHECODES"]
        kind_label = {
            "CONTINUOUS_VARIABLE": "Continuous variables",
            "CV_ENDPOINTS": "Cardiovascular endpoints",
            "SELF_REPORTED": "Self reported",
            "PHECODES": "Phecodes",
        }
        plot_df = plot_df[plot_df["analysis_type"].isin(kind_order)].copy()
        if plot_df.empty:
            fig, ax = plt.subplots()
            ax.text(0.5, 0.5, "No rows with recognized analysis_type.", ha="center", va="center")
            ax.set_axis_off()
            return fig

        x_pos_of = {k: i for i, k in enumerate(kind_order)}
        plot_df["x"] = plot_df["analysis_type"].map(x_pos_of)

        # Compute y and jitter per row; reset index so group indices match arrays
        plot_df = plot_df.reset_index(drop=True)
        if use_log:
            plot_df["_y"] = -np.log10(plot_df[metric].astype(float))
            y_label = f"−log10({metric})"
            thresh_y = -math.log10(0.05)
        else:
            plot_df["_y"] = plot_df[metric].astype(float)
            y_label = metric
            thresh_y = 0.05

        rng = np.random.default_rng(0)
        plot_df["_xj"] = plot_df["x"].astype(float) + rng.normal(0, 0.045, size=len(plot_df))

        # Plot per gene using per-row columns (no external array indexing)
        fig, ax = plt.subplots(figsize=(7.2, 3.8))
        for gene_name, gdf in plot_df.groupby("gene", sort=True):
            ax.scatter(gdf["_xj"], gdf["_y"], s=20, alpha=0.8, label=gene_name)

        ax.axhline(thresh_y, linestyle=":", linewidth=1)
        ax.set_xticks(range(len(x_pos_of)))
        ax.set_xticklabels([kind_label[k] for k in kind_order])
        ax.set_ylabel(y_label)
        ax.set_title(f"{y_label} by analysis type (top {limit}/gene/group)")

        if input.show_legend():
            ax.legend(title="Gene", loc="upper right", frameon=True, fontsize="small")

        fig.tight_layout()
        return fig

    def get_sub_df(kind: str) -> pd.DataFrame:
        df = df_results.get()
        if df is None or df.empty or "analysis_type" not in df.columns:
            return tidy_for_table(pd.DataFrame())
        sub = df[df["analysis_type"] == kind].copy()
        sort_cols = [c for c in ["q", "p"] if c in sub.columns]
        if sort_cols:
            sub = sub.sort_values(sort_cols, ascending=[True] * len(sort_cols))
        if "gene" not in sub.columns:
            sub["gene"] = "unknown"
        N = max(1, int(input.limit()))
        sub = sub.groupby("gene", group_keys=False).head(N)
        return tidy_for_table(sub)

    @output
    @render.table
    def tbl_continuous():
        return get_sub_df("CONTINUOUS_VARIABLE")

    @output
    @render.table
    def tbl_cv():
        return get_sub_df("CV_ENDPOINTS")

    @output
    @render.table
    def tbl_self():
        return get_sub_df("SELF_REPORTED")

    @output
    @render.table
    def tbl_phecode():
        return get_sub_df("PHECODES")

    @output
    @render.text
    def log_text():
        return logs.get() or "…"

app = App(app_ui, server)
