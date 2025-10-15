import sys
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import shiny.ui as ui
from shiny import App, reactive, render

from backend import (
    fetch_gene_results,
    get_continuous_labels,
    get_cv_labels,
    get_phecode_labels,
    get_self_reported_labels,
    get_single_gene_df,
    outcome_catalog,
    parse_gene_list,
    resolve_gene,
    tidy_table,
    top_n_per_gene,
)
from plotting import (
    add_jitter,
    apply_scale,
    bar_plot,
    bubble_plot,
    heatmap_plot,
    prepare_plot_df,
    volcano_plot,
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
    nan_color = "#FF00FF"

    ###########################################################################
    # Utility / helper functions
    ###########################################################################
    def _tbl(kind: str, filters: list):
        """Tables: Get top N genes per endpoint"""
        # Get full dataframe
        df = df_results.get()

        # If we want to copy and paste from a list rather than use presets
        if filters:
            if not input.preset():
                filters = filters.split("\n")

        # No data loaded or empty data frame of results
        if df is None or df.empty:
            return tidy_table(
                pd.DataFrame(),
                metric=input.metric(),
                threshold=float(input.threshold()),
                filters=None,
            )

        # If data loaded and non-empty, we limit to the desired number of results
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
        """Filter either by p or q value (FDR)"""
        if str(input.metric()) in df.columns:
            df = df[df[str(input.metric())] < float(input.threshold())]
        return df

    def _safe_input(input, name):
        """Safe get conditional input panels which haven't been added to the UI"""
        try:
            return getattr(input, name)()
        except AttributeError:
            print(f"Input {name} not yet created on UI", file=sys.stderr)
            return None

    def _make_selectize(id_, label, choices, preset=True):
        if preset:
            return ui.input_selectize(
                id_,
                label,
                choices=choices,
                multiple=True,
                selected=sel_keep.get(id_, []),
                options={
                    "create": True,
                    "persist": False,
                    "placeholder": "Type to search…",
                },
            )
        else:
            print("adding text area")
            return ui.input_text_area(id_, label)

    ###########################################################################
    # Effects and events
    ###########################################################################
    sel_keep = {  # should not be reactive, otherwise we have a classical reactive infinite loop
        "filter_cont": [],
        "filter_cv": [],
        "filter_self": [],
        "filter_phe": [],
    }

    @reactive.Effect
    @reactive.event(input.nan_color)
    def _nan_color():
        nan_color = input["nan_color"]()
        print(f"nan_color: {nan_color}")

    @reactive.Effect
    @reactive.event(input.filter_cont)
    def _save_cont():
        sel_keep["filter_cont"] = input.filter_cont() or []

    @reactive.Effect
    @reactive.event(input.filter_cv)
    def _save_cv():
        sel_keep["filter_cv"] = input.filter_cv() or []

    @reactive.Effect
    @reactive.event(input.filter_self)
    def _save_self():
        sel_keep["filter_self"] = input.filter_self() or []

    @reactive.Effect
    @reactive.event(input.filter_phe)
    def _save_phe():
        sel_keep["filter_phe"] = input.filter_phe() or []

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
    @render.download(filename="results.csv")
    def download_csv():
        """Download helper"""
        df = df_results.get()
        if df is None or df.empty:
            dl_msg.set("No data to download.")
            return
        # TODO: Yield all data currently, perhaps want filtering here too before download?
        _log("User downloaded CSV file")
        yield df.to_csv(index=False).encode("utf-8")

    @output
    @render.text
    def download_msg():
        """Status message for download"""
        return dl_msg.get()

    @reactive.Calc
    def prepared_df():
        """Cache preparation of df"""
        # get dataframe of results
        df = df_results.get()
        # global filter for p or q value
        df = _filter_for_p_or_q_value(df)
        return df

    @output
    @render.plot
    def plot_out():
        # Get cache df
        df = prepared_df()

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

        filters = {
            "filter_cont": (
                _safe_input(input, "filter_cont") if input.use_cont() else None
            ),
            "filter_cv": _safe_input(input, "filter_cv") if input.use_cv() else None,
            "filter_self": (
                _safe_input(input, "filter_self") if input.use_self() else None
            ),
            "filter_phe": _safe_input(input, "filter_phe") if input.use_phe() else None,
        }

        # if we do not with to use select from database
        if not input.preset():
            filters["filter_cont"] = (
                input.filter_cont().split("\n") if input.use_cont() else None
            )
            filters["filter_cv"] = (
                input.filter_cv().split("\n") if input.use_cv() else None
            )
            filters["filter_self"] = (
                input.filter_self().split("\n") if input.use_self() else None
            )
            filters["filter_phe"] = (
                input.filter_phe().split("\n") if input.use_phe() else None
            )

        d = prepare_plot_df(df, metric, limit, outcome_catalog(), filters)

        print("prepared plot:")
        print(d)

        if d.empty or input.plot_type() is None:
            fig, ax = plt.subplots()
            ax.text(0.5, 0.5, "No rows to plot.", ha="center", va="center")
            ax.set_axis_off()
            return fig

        # prepare plot
        d = apply_scale(d, use_log, metric, eps=1e-300)
        d = add_jitter(d, seed=0, sd=0.045)

        # plot with given type
        if str(input.plot_type()) == "Volcano plot":
            return volcano_plot(
                d,
                int(input.limit()),
                bool(input.neglog10()),
                bool(input.show_legend()),
                str(input.metric()),
            )

        if str(input.plot_type()) == "Bar plot":
            return bar_plot(d, input.show_legend())

        if str(input.plot_type()) == "Heatmap":
            print("nan_color:")
            print(nan_color)
            return heatmap_plot(
                d,
                input.metric(),
                bool(input.neglog10()),
                float(input.threshold()),
                str(nan_color),
            )
        if str(input.plot_type()) == "Bubble plot":
            gene = input.single_gene().upper()
            category = input.single_gene_category()
            return bubble_plot(
                get_single_gene_df(d, gene, category),
                category,
                gene,
                input.metric(),
                input.show_legend(),
            )

    @output
    @render.data_frame
    def tbl_continuous():
        """CONTINUOUS VARIABLES"""
        if input.use_cont():
            return _tbl("CONTINUOUS_VARIABLE", _safe_input(input, "filter_cont"))
        return _tbl("CONTINUOUS_VARIABLE", None)

    @output
    @render.data_frame
    def tbl_cv():
        """CARDIOVASCULAR VARIABLES"""
        if input.use_cv():
            return _tbl("CV_ENDPOINTS", _safe_input(input, "filter_cv"))
        return _tbl("CV_ENDPOINTS", None)

    @output
    @render.data_frame
    def tbl_self():
        """SELF REPORTED"""
        if input.use_self():
            return _tbl("SELF_REPORTED", _safe_input(input, "filter_self"))
        return _tbl("SELF_REPORTED", None)

    @output
    @render.data_frame
    def tbl_phecode():
        """PHECODES"""
        if input.use_phe():
            return _tbl("PHECODES", _safe_input(input, "filter_phe"))
        return _tbl("PHECODES", None)

    @output
    @render.text
    def log_text():
        """ " Log message"""
        return logs.get() or "…"

    ############################################################################
    # DYNAMIC UI CONTENT
    ############################################################################
    @output
    @render.ui
    def filters_nav():
        panels = []
        if input.use_cont():
            panels.append(
                ui.nav_panel(
                    "Continuous variables",
                    _make_selectize(
                        "filter_cont",
                        "Select continuous variables",
                        get_continuous_labels(),
                        input.preset(),
                    ),
                )
            )
        if input.use_cv():
            panels.append(
                ui.nav_panel(
                    "CV endpoints",
                    _make_selectize(
                        "filter_cv",
                        "Select CV endpoints",
                        get_cv_labels(),
                        input.preset(),
                    ),
                )
            )
        if input.use_self():
            panels.append(
                ui.nav_panel(
                    "Self reported",
                    _make_selectize(
                        "filter_self",
                        "Select self reported",
                        get_self_reported_labels(),
                        input.preset(),
                    ),
                )
            )
        if input.use_phe():
            panels.append(
                ui.nav_panel(
                    "Phecodes",
                    _make_selectize(
                        "filter_phe",
                        "Select phecodes",
                        get_phecode_labels(),
                        input.preset(),
                    ),
                )
            )

        if not panels:
            return ui.div(
                {"class": "text-muted"},
                "Enable filter categories with the checkboxes above.",
            )

        return ui.navset_pill(*panels)


###############################################################################
# START APPLICATION
###############################################################################
STATIC_DIR = Path(__file__).parent / "www"
app = App(app_ui, server, static_assets=STATIC_DIR)
