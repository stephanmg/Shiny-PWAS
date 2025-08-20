import math

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from shiny import App, reactive, render

from backend import (
    fetch_gene_results,
    get_single_gene_df,
    outcome_catalog,
    parse_gene_list,
    resolve_gene,
    tidy_table,
    top_n_per_gene,
)
from plotting import KIND_LABEL, KIND_ORDER, add_jitter, apply_scale, prepare_plot_df
from ui import make_ui

###############################################################################
# UI
###############################################################################
app_ui = make_ui()


def four_heatmaps(df, metric="p", use_log=True, pthresh=None):
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
    fig, axes = plt.subplots(2, 2, figsize=(12, 8))
    axes = axes.ravel()

    vmin, vmax = d["_val"].min(), d["_val"].max()

    for ax, kind in zip(axes, KIND_ORDER):
        sub = d[d["analysis_type"] == kind]
        if sub.empty:
            ax.text(
                0.5, 0.5, f"No data for\n{KIND_LABEL[kind]}", ha="center", va="center"
            )
            ax.set_axis_off()
            continue

        # pivot: rows=Description, cols=genes
        mat = sub.pivot(index="Description", columns="gene", values="_val").fillna(0)

        # set bad color to magenta (missing)
        cmap = plt.cm.viridis.copy()
        cmap.set_bad(color="magenta")

        im = ax.imshow(
            mat.values,
            aspect="auto",
            interpolation="nearest",
            vmin=vmin,
            vmax=vmax,
            cmap="viridis",
        )

        ax.set_title(KIND_LABEL[kind])
        ax.set_xticks(range(len(mat.columns)))
        ax.set_xticklabels(mat.columns, rotation=45, ha="right")
        ax.set_yticks(range(len(mat.index)))
        ax.set_yticklabels(mat.index, fontsize=6)

        cbar = fig.colorbar(im, ax=ax, shrink=0.7)
        cbar.set_label(f"-log10({metric})" if use_log else metric)

    return fig


###############################################################################
# Server / business logic
###############################################################################
def server(input, output, session):
    """Server"""
    logs = reactive.Value("")
    df_results = reactive.Value(pd.DataFrame())
    dl_msg = reactive.Value(
        "Click 'Download CSV' to save the combined table with all rows."
    )

    ###########################################################################
    # Utility / helper functions
    ###########################################################################
    def _tbl(kind: str, filters: list):
        """Tables: Get top N genes per endpoint"""
        # Get full dataframe
        df = df_results.get()

        if df is None or df.empty:
            return tidy_table(
                pd.DataFrame(),
                metric=input.metric(),
                threshold=float(input.threshold()),
                filters=None,
            )
        sub = top_n_per_gene(df, kind, int(input.limit()))
        return tidy_table(
            sub,
            metric=input.metric(),
            threshold=float(input.threshold()),
            filters=filters,
        )

    def _log(msg):
        """Log message to logging window"""
        logs.set((logs.get() or "") + msg + "\n")

    def _filter_for_p_or_q_value(df: pd.DataFrame):
        if str(input.metric()) in df.columns:
            df = df[df[str(input.metric())] < float(input.threshold())]
        return df

    def _filter_for_phenotypes_by_category(df: pd.DataFrame):
        mask = pd.Series(True, index=df.index)
        mask &= df["Description"].isin(input.filter_cont())
        return mask

    ###########################################################################
    # Effects and events
    ###########################################################################
    @reactive.Effect
    @reactive.event(input.btn_phenos, input.genes)
    def load_phenos():
        """Load phenotypes for selected gene from database ExPheWAS"""
        genes = parse_gene_list(input.genes() or "")
        subset = input.subset()
        if not genes:
            _log("Enter at least one gene (symbol or ENSG).")
            df_results.set(pd.DataFrame())
            return
        frames = []
        for g in genes:
            ensg, sym = resolve_gene(g)
            if not ensg:
                _log(f"! Could not find gene with symbol '{g}' — skipping.")
                continue
            try:
                df = fetch_gene_results(ensg, subset)
                if not df.empty:
                    df = df.copy()
                    df["gene"] = sym
                    frames.append(df)
                    _log(f"{sym}: {len(df)} rows.")
                else:
                    _log(f"No results found for gene symbol {sym}.")
            except Exception as e:
                _log(f"{sym}: ERROR {e}")
        if frames:
            all_df = pd.concat(frames, ignore_index=True)
            df_results.set(all_df)
            _log(f"Total rows combined: {len(all_df)}")
        else:
            df_results.set(pd.DataFrame())
            _log("No data loaded.")

    @output
    @render.download(filename="results.csv")
    def download_csv():
        """Download helper"""
        df = df_results.get()
        if df is None or df.empty:
            dl_msg.set("No data to download.")
            return
        tidy = tidy_table(df.copy())
        yield tidy.to_csv(index=False).encode("utf-8")
        _log("User downloaded CSV file to {filename}")

    @reactive.Effect
    @reactive.event(input.btn_plot)
    def log_plot():
        _log("Plot requested.")

    @reactive.Effect
    @reactive.event(input.threshold)
    def log_threshold():
        _log(f"Threshold requested: {input.threshold()}")

    ###########################################################################
    # Outputs
    ###########################################################################
    @output
    @render.text
    def download_msg():
        return dl_msg.get()

    @output
    @render.plot
    def plot_out():
        df = df_results.get()
        df = _filter_for_p_or_q_value(df)

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
        if (
            len(
                input.filter_cont()
                + input.filter_cv()
                + input.filter_self()
                + input.filter_phe()
            )
            == 0
        ):
            d = prepare_plot_df(df, metric, limit, outcome_catalog(), None)
        else:
            filters = {
                "filter_cont": input.filter_cont(),
                "filter_cv": input.filter_cv(),
                "filter_self": input.filter_self(),
                "filter_phe": input.filter_phe(),
            }
            d = prepare_plot_df(df, metric, limit, outcome_catalog(), filters)

        if d.empty or input.plot_type() is None:
            fig, ax = plt.subplots()
            ax.text(0.5, 0.5, "No rows to plot.", ha="center", va="center")
            ax.set_axis_off()
            return fig

        # prepare plot
        d = apply_scale(d, use_log, metric, eps=1e-300)
        d = add_jitter(d, seed=0, sd=0.045)

        if str(input.plot_type()) == "Volcano plot":
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
                ax.legend(
                    title="Gene", loc="upper right", frameon=True, fontsize="small"
                )
            fig.tight_layout()
            return fig

        if str(input.plot_type()) == "Bar plot":
            # Ensure categories are consistent
            order = ["CONTINUOUS_VARIABLE", "CV_ENDPOINTS", "SELF_REPORTED", "PHECODES"]
            labels = {
                "CONTINUOUS_VARIABLE": "Continuous",
                "CV_ENDPOINTS": "CV endpoints",
                "SELF_REPORTED": "Self reported",
                "PHECODES": "Phecodes",
            }

            # Count rows per (gene, analysis_type)
            counts = (
                d.groupby(["gene", "analysis_type"]).size().reset_index(name="count")
            )

            genes = sorted(d["gene"].unique())
            x = np.arange(len(genes))
            width = 0.2

            fig, ax = plt.subplots(figsize=(8, 4.5))

            for i, kind in enumerate(order):
                subset = counts[counts["analysis_type"] == kind]
                # Align with all genes (fill missing with 0)
                y = [subset.loc[subset["gene"] == g, "count"].sum() for g in genes]
                ax.bar(x + i * width, y, width, label=labels[kind])

            ax.set_xticks(x + width * (len(order) - 1) / 2)
            ax.set_xticklabels(genes, rotation=45, ha="right")
            ax.set_ylabel("Number of outcomes")
            ax.set_title("Counts per gene by analysis type")
            ax.legend(title="Category")
            fig.tight_layout()
            return fig

        if str(input.plot_type()) == "Heatmap":
            return four_heatmaps(
                d, input.metric(), bool(input.neglog10()), float(input.threshold())
            )
        if str(input.plot_type()) == "Bubble plot":
            gene = input.single_gene()
            category = input.single_gene_category()
            d = get_single_gene_df(d, gene, category)
            if d.empty:
                fig, ax = plt.subplots()
                ax.text(0.5, 0.5, "No rows for single gene.", ha="center", va="center")
                ax.set_axis_off()
                return fig

            sizes = d["n_cases"] / d["n_cases"].max() * 800  # adjust scaling factor

            metric = input.metric()
            if metric not in d.columns:
                metric = "p"  # fallback
            vals = d[metric].replace(0, 1e-300)  # avoid log(0)

            # Use -log10 for better spread
            colors = -np.log10(vals)

            fig, ax = plt.subplots(figsize=(10, 5))
            sc = ax.scatter(
                d["n_controls"],
                d["Description"],
                s=sizes,
                c=colors,
                cmap="viridis_r",
                alpha=0.8,
                edgecolor="k",
            )

            ax.set_xlabel("Number of controls")
            ax.set_ylabel("Phenotype (Description)")
            ax.set_title(f"Bubble plot for {gene} (colored by −log10({metric}))")
            ax.grid(True, linestyle="--", alpha=0.5)
            ax.tick_params(axis="y", labelsize=8)

            # Colorbar for metric
            cbar = fig.colorbar(sc, ax=ax, shrink=0.8)
            cbar.set_label(f"−log10({metric})")

            # Legend for bubble sizes
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

    @output
    @render.table
    def tbl_continuous():
        return _tbl("CONTINUOUS_VARIABLE", input.filter_cont())

    @output
    @render.table
    def tbl_cv():
        return _tbl("CV_ENDPOINTS", input.filter_cv())

    @output
    @render.table
    def tbl_self():
        return _tbl("SELF_REPORTED", input.filter_self())

    @output
    @render.table
    def tbl_phecode():
        return _tbl("PHECODES", input.filter_phe())

    @output
    @render.text
    def log_text():
        return logs.get() or "…"


###############################################################################
# Start application
###############################################################################
app = App(app_ui, server)
