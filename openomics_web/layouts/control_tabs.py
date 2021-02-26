import dash_core_components as dcc
import dash_html_components as html

from openomics import MicroRNA, LncRNA, MessengerRNA, Protein


def control_tabs():
    return html.Div(
        id="circos-control-tabs",
        className="control-tabs",
        children=[
            dcc.Tabs(
                id="circos-tabs",
                value="what-is",
                children=[
                    dcc.Tab(
                        label="Welcome",
                        value="what-is",
                        children=html.Div(
                            className="control-tab",
                            children=[
                                html.H5(
                                    className="what-is",
                                    children=
                                    "OpenOmics: The Multi-Omics Explorer",
                                ),
                                html.
                                P("OpenOmics provides a bioinformatics API and web-app "
                                  "platform integrate, analyze, and visualize the "
                                  "multi-omics and clinical data. You can interact with the "
                                  "data input tabs to integrate the various different data types "
                                  "including bio-databases, multi-omics expression, genomics, "
                                  "and clinical data."),
                                html.Div([
                                    "Cite us at: ",
                                    html.A(
                                        "OpenOmics paper",
                                        href=
                                        "http://www.doi.org/10.1101/gr.092759.109",
                                    ),
                                    ".",
                                ]),
                                html.Div([
                                    "For a look into the OpenOmics API and usage documentation, please visit the "
                                    "original repository ",
                                    html.A(
                                        "OpenOmics GitHub",
                                        href=
                                        "https://github.com/JonnyTran/OpenOmics.git",
                                    ),
                                    ".",
                                ]),
                                html.Br(),
                                html.Button(
                                    id="start-button",
                                    style={"align": "center"},
                                    n_clicks=0,
                                    children="Get Started",
                                ),
                            ],
                        ),
                    ),
                    dcc.Tab(
                        label="Import Clinical Data",
                        value="clinical-upload",
                        className="control-tab",
                        children=[
                            html.Div([
                                "Cohort name: ",
                            ]),
                            dcc.Input(
                                id="clinical-cohort",
                                style={"width": "100%"},
                                placeholder="Leave empty to use filename",
                                type="text",
                            ),
                            html.Br(),
                            html.Div(["Select clinical data type"]),
                            dcc.RadioItems(
                                id="clinical-data-type",
                                options=[
                                    {
                                        "label": "Patients data",
                                        "value": "Patients data",
                                    },
                                    {
                                        "label": "Samples data",
                                        "value": "Samples data"
                                    },
                                    {
                                        "label": "Drug responses",
                                        "value": "Drug responses",
                                    },
                                ],
                                value="Patients data",
                            ),
                            html.Br(),
                            html.Div([
                                "Import Clinical Data",
                            ]),
                            html.Div(children=[
                                dcc.Upload(
                                    id="upload-clinical",
                                    children=html.Div([
                                        "Drag and Drop or ",
                                        html.A("Select Files"),
                                    ]),
                                    style={
                                        "width": "100%",
                                        "height": "60px",
                                        "lineHeight": "60px",
                                        "borderWidth": "1px",
                                        "borderStyle": "dashed",
                                        "borderRadius": "5px",
                                        "textAlign": "center",
                                        "margin": "10px",
                                    },
                                    multiple=False,
                                )
                            ]),
                            html.Br(),
                            html.Div(
                                id="clinical-column-select",
                                children=[
                                    html.Div([
                                        "Select the patient id column to index by:"
                                    ]),
                                    dcc.Dropdown(
                                        id="clinical-patient-col-name",
                                        disabled=True),
                                    html.Div([
                                        "Select the column prefixes to import:"
                                    ]),
                                    dcc.Dropdown(
                                        id="clinical-data-columns-select",
                                        disabled=True,
                                    ),
                                ],
                            ),
                            html.Br(),
                            html.Button(
                                id="clinical-submit-button",
                                n_clicks=0,
                                children="Submit",
                            ),
                        ],
                    ),
                    dcc.Tab(
                        label="Import Omics Data",
                        value="data-upload",
                        className="control-tab",
                        children=[
                            html.Div([
                                "Cohort name: ",
                            ]),
                            dcc.Input(
                                id="data-table-cohort",
                                style={"width": "100%"},
                                placeholder="Leave empty to use filename",
                                type="text",
                            ),
                            html.Br(),
                            html.Br(),
                            html.Div([
                                "Select data type:",
                            ]),
                            dcc.RadioItems(
                                id="data-table-type",
                                options=[
                                    {
                                        "label": "Protein Expression",
                                        "value": Protein.name(),
                                    },
                                    {
                                        "label": "miRNA Expression",
                                        "value": MicroRNA.name(),
                                    },
                                    {
                                        "label": "lncRNA Expression",
                                        "value": LncRNA.name(),
                                    },
                                    {
                                        "label": "mRNA Expression",
                                        "value": MessengerRNA.name(),
                                    },
                                ],
                                value=MicroRNA.name(),
                            ),
                            html.Br(),
                            html.Div([
                                "Import a data file",
                            ]),
                            html.Div(children=[
                                dcc.Upload(
                                    id="upload-data-table",
                                    children=html.Div([
                                        "Drag and Drop or ",
                                        html.A("Select Files"),
                                    ]),
                                    style={
                                        "width": "100%",
                                        "height": "60px",
                                        "lineHeight": "60px",
                                        "borderWidth": "1px",
                                        "borderStyle": "dashed",
                                        "borderRadius": "5px",
                                        "textAlign": "center",
                                        "margin": "10px",
                                    },
                                    multiple=True,
                                )
                            ]),
                            html.Br(),
                            html.Div(
                                id="data-table-column-select",
                                children=[
                                    html.Div([
                                        "Select the gene id/name column to index by:"
                                    ]),
                                    dcc.Dropdown(
                                        id="data-table-genes-col-name",
                                        disabled=True),
                                    html.Div([
                                        "Select the column prefixes to import:"
                                    ]),
                                    dcc.Dropdown(
                                        id="data-table-columns-select",
                                        disabled=True,
                                    ),
                                ],
                            ),
                            html.Br(),
                            html.Div(["Each gene is indexed by:"]),
                            dcc.RadioItems(
                                id="data-table-genes_index",
                                options=[
                                    {
                                        "label": "Gene symbol",
                                        "value": "gene_name"
                                    },
                                    {
                                        "label": "Gene index",
                                        "value": "gene_index"
                                    },
                                    {
                                        "label": "Transcript name",
                                        "value": "transcript_name",
                                    },
                                    {
                                        "label": "Transcript index",
                                        "value": "transcript_index",
                                    },
                                ],
                                value="gene_name",
                            ),
                            html.Br(),
                            html.Div(["Is the table transposed?"]),
                            dcc.RadioItems(
                                id="data-table-transpose",
                                options=[
                                    {
                                        "label": "Sample name columns",
                                        "value": "True"
                                    },
                                    {
                                        "label": "Gene name in columns",
                                        "value": "False"
                                    },
                                ],
                                value="True",
                            ),
                            html.Br(),
                            html.Button(
                                id="upload-data-table-submit",
                                n_clicks=0,
                                children="Submit",
                            ),
                        ],
                    ),
                    dcc.Tab(
                        label="Annotate database",
                        value="annotatte-db",
                        className="control-tab",
                        children=[
                            html.Div("To be implemented"),
                        ],
                    ),
                ],
            )
        ],
    )
