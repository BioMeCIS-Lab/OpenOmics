import dash_core_components as dcc
import dash_html_components as html

from openomics_web.layouts.control_tabs import control_tabs


def app_main():
    return html.Div(
        id="circos-body",
        className="app-body",
        children=[
            loading(),
            control_tabs(),
            html.Div(id="output-data-upload", className="circos-size"),
            html.Div(id="output-clinical-upload", className="circos-size"),
        ],
    )


def loading():
    return html.Div(
        id="loading",
        children=[
            dcc.Loading(
                className="loading",
                children=html.Div(id="circos-hold",
                                  # children=[empty]
                                  ),
            )
        ],
    )
