import pandas as pd
import plotly.graph_objects as go
from shiny import App, reactive, render
from shinywidgets import render_plotly

from backend import (
    fetch_gene_results,
    outcome_catalog,
    parse_gene_list,
    resolve_gene,
    tidy_table,
    top_n_per_gene,
)
from plotting import (
    add_jitter,
    apply_scale,
    plot_bar,
    plot_bubble,
    plot_heatmap,
    plot_volcano,
    prepare_plot_df,
)
from ui import make_ui

###############################################################################
# UI
###############################################################################
app_ui = make_ui()


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

    @render_plotly
    def plot_out():
        df = df_results.get()
        df = _filter_for_p_or_q_value(df)

        if df is None or df.empty:
            return go.Figure().add_annotation(
                text="No data yet. Click 'Load phenotypes'.",
                showarrow=False,
                x=0.5,
                y=0.5,
            )

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

        if d is None or d.empty:
            fig = go.Figure()
            fig.add_annotation(text="No rows to plot.", x=0.5, y=0.5, showarrow=False)
            fig.update_xaxes(visible=False)
            fig.update_yaxes(visible=False)
            fig.update_layout(margin=dict(l=10, r=10, t=30, b=10))
            return fig

        # prepare plot further
        d = apply_scale(d, use_log, metric, eps=1e-300)
        d = add_jitter(d, seed=0, sd=0.045)

        # decide on visualization
        if str(input.plot_type()) == "Volcano plot":
            return plot_volcano(d, use_log, metric, limit, input.show_legend())

        if str(input.plot_type()) == "Bar plot":
            return plot_bar(d)

        if str(input.plot_type()) == "Heatmap":
            return plot_heatmap(
                d, input.metric(), bool(input.neglog10()), float(input.threshold())
            )

        if str(input.plot_type()) == "Bubble plot":
            return plot_bubble(
                d, input.single_gene(), input.single_gene_category(), input.metric()
            )

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
