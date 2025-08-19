from shiny import ui


def make_ui():
    return ui.page_sidebar(
        #######################################################################
        # SIDEBAR (LEFT)
        #######################################################################
        ui.sidebar(
            ui.input_text_area(
                "genes",
                "Gene names (comma or newline separated)",
                value="METTL2A\nPCSK9",
                width="100%",
                height="120px",
            ),
            ui.card(
                ui.card_header("Filter phenotypes by category"),
                ui.navset_pill(
                    ui.nav_panel(
                        "Continuous variables",
                        ui.input_selectize(
                            "filter_cont",
                            "Select continuous variables",
                            choices=[],
                            multiple=True,
                            options={
                                "create": True,
                                "persist": False,
                                "placeholder": "Type to search…",
                            },
                        ),
                    ),
                    ui.nav_panel(
                        "Cardiovascular endpoints",
                        ui.input_selectize(
                            "filter_cv",
                            "Select cardiovascular endpoints",
                            choices=[],
                            multiple=True,
                            options={
                                "create": True,
                                "persist": False,
                                "placeholder": "Type to search…",
                            },
                        ),
                    ),
                    ui.nav_panel(
                        "Self reported",
                        ui.input_selectize(
                            "filter_self",
                            "Select self reported phenotypes",
                            choices=[],
                            multiple=True,
                            options={
                                "create": True,
                                "persist": False,
                                "placeholder": "Type to search…",
                            },
                        ),
                    ),
                    ui.nav_panel(
                        "Phecodes",
                        ui.input_selectize(
                            "filter_phe",
                            "Select phecodes",
                            choices=[],
                            multiple=True,
                            options={
                                "create": True,
                                "persist": False,
                                "placeholder": "Type to search…",
                            },
                        ),
                    ),
                ),
            ),
            ui.input_radio_buttons(
                "subset",
                "Subset selection",
                choices={
                    "BOTH": "Both",
                    "MALE_ONLY": "Male only",
                    "FEMALE_ONLY": "Female only",
                },
                selected="BOTH",
            ),
            ui.input_numeric(
                "limit", "Max rows per gene per group", value=1e6, min=1, step=1
            ),
            ui.input_action_button("btn_phenos", "Load phenotypes"),
            ui.input_select(
                "metric",
                "Value to plot",
                choices={"p": "p-value", "q": "q-value (FDR)"},
                selected="p",
            ),
            ui.input_select(
                "threshold",
                "Threshold",
                choices={"0.05": "0.05", "0.01": "0.01", "0.001": "0.001"},
                selected="0.05",
            ),
            ui.input_checkbox("neglog10", "Use −log10 scale", value=True),
            ui.input_checkbox("show_legend", "Show legend", value=True),
            ui.input_radio_buttons(
                "plot_type",
                "Visualization type",
                choices={
                    "Volcano plot": "Volcano plot",
                    "Bar plot": "Bar plot",
                    "Heatmap": "Heatmap",
                    "Bubble plot": "Bubble plot",
                },
                selected="Volcano plot",
            ),
            ui.panel_conditional(
                "input.plot_type == 'Bubble plot'",
                ui.input_text_area("single_gene", "Single gene", value="PCSK9"),
                ui.input_selectize(
                    "single_gene_category",
                    "Select category",
                    choices=["CV_ENDPOINTS", "SELF_REPORTED", "PHECODES"],
                ),
            ),
            ui.input_action_button("btn_plot", "Load plot"),
            width=340,
        ),
        #######################################################################
        # MAIN PANEL (RIGHT)
        #######################################################################
        ui.tags.style(".content-wrapper {max-width: 1100px; margin: 0 auto;}"),
        ui.div(
            {"class": "content-wrapper"},
            ui.card(
                ui.card_header("Plot of selected metric by analysis type"),
                ui.output_plot("plot_out", height="380px"),
            ),
            ui.card(
                ui.card_header("Phenotypes (top N genes)"),
                ui.navset_pill(
                    ui.nav_panel(
                        "Continuous variables (CVs)", ui.output_table("tbl_continuous")
                    ),
                    ui.nav_panel("Cardiovascular endpoints", ui.output_table("tbl_cv")),
                    ui.nav_panel("Self reported", ui.output_table("tbl_self")),
                    ui.nav_panel("Phecodes", ui.output_table("tbl_phecode")),
                ),
            ),
            ui.card(
                ui.card_header("Download combined table (all rows)"),
                ui.input_action_button("btn_download_csv", "Download CSV"),
                ui.output_text_verbatim("download_msg"),
            ),
            ui.card(
                ui.card_header("Log"),
                ui.output_text_verbatim("log_text"),
            ),
        ),
        #######################################################################
        # APP TITLE
        #######################################################################
        title="Shiny-PWAS: multi-gene phenotypes view",
    )
