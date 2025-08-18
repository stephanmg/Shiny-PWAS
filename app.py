# app.py
import math

import matplotlib.pyplot as plt
import pandas as pd
from shiny import App, reactive, render

from backend import (
    fetch_gene_results,
    parse_gene_list,
    resolve_gene,
    tidy_table,
    top_n_per_gene,
)
from plotting import KIND_LABEL, KIND_ORDER, add_jitter, apply_scale, prepare_plot_df
from ui import make_ui

app_ui = make_ui()


def server(input, output, session):
    """Server"""
    logs = reactive.Value("")
    df_results = reactive.Value(pd.DataFrame())
    dl_msg = reactive.Value(
        "Click 'Download CSV' to save the combined table with all rows."
    )

    def log(msg):
        """Log message to logging window"""
        logs.set((logs.get() or "") + msg + "\n")

    @reactive.Effect
    @reactive.event(input.btn_phenos)
    def load_phenos():
        """Load phenotypes for selected gene from database ExPheWAS"""
        genes = parse_gene_list(input.genes() or "")
        subset = input.subset()
        if not genes:
            log("Enter at least one gene (symbol or ENSG).")
            df_results.set(pd.DataFrame())
            return
        frames = []
        for g in genes:
            ensg, sym = resolve_gene(g)
            if not ensg:
                log(f"! Could not resolve '{g}' — skipping.")
                continue
            try:
                df = fetch_gene_results(ensg, subset)
                if not df.empty:
                    df = df.copy()
                    df["gene"] = sym
                    frames.append(df)
                    log(f"{sym}: {len(df)} rows.")
                else:
                    log(f"{sym}: 0 rows.")
            except Exception as e:
                log(f"{sym}: ERROR {e}")
        if frames:
            all_df = pd.concat(frames, ignore_index=True)
            df_results.set(all_df)
            log(f"Total rows combined: {len(all_df)}")
        else:
            df_results.set(pd.DataFrame())
            log("No data loaded.")

    @reactive.Effect
    @reactive.event(input.btn_download_csv)
    def dl():
        """Download helper"""
        df = df_results.get()
        if df is None or df.empty:
            dl_msg.set("No data to download.")
            return
        tidy = tidy_table(df.copy())
        path = "/mnt/data/multigene_results_all.csv"
        tidy.to_csv(path, index=False)
        dl_msg.set(f"Saved CSV to {path}")

    @output
    @render.text
    def download_msg():
        return dl_msg.get()

    @reactive.Effect
    @reactive.event(input.btn_plot)
    def log_plot():
        log("Plot requested.")

    @output
    @render.plot
    def plot_out():
        df = df_results.get()
        if df is None or df.empty:
            fig, ax = plt.subplots()
            ax.text(
                0.5,
                0.5,
                "No data yet. Click 'Load phenotypes'.",
                ha="center",
                va="center",
            )
            ax.set_axis_off()
            return fig
        metric = input.metric()
        use_log = bool(input.neglog10())
        limit = max(1, int(input.limit()))
        d = prepare_plot_df(df, metric, limit)
        if d.empty:
            fig, ax = plt.subplots()
            ax.text(0.5, 0.5, "No rows to plot.", ha="center", va="center")
            ax.set_axis_off()
            return fig
        d = apply_scale(d, use_log, metric, eps=1e-300)
        d = add_jitter(d, seed=0, sd=0.045)

        # draw
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
        if input.show_legend():
            ax.legend(title="Gene", loc="upper right", frameon=True, fontsize="small")
        fig.tight_layout()
        return fig

    # tables: top N per gene, per endpoint
    def _tbl(kind: str):
        df = df_results.get()
        if df is None or df.empty:
            return tidy_table(pd.DataFrame())
        sub = top_n_per_gene(df, kind, int(input.limit()))
        return tidy_table(sub)

    @output
    @render.table
    def tbl_continuous():
        return _tbl("CONTINUOUS_VARIABLE")

    @output
    @render.table
    def tbl_cv():
        return _tbl("CV_ENDPOINTS")

    @output
    @render.table
    def tbl_self():
        return _tbl("SELF_REPORTED")

    @output
    @render.table
    def tbl_phecode():
        return _tbl("PHECODES")

    @output
    @render.text
    def log_text():
        return logs.get() or "…"


app = App(app_ui, server)
