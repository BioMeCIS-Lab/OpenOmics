import dash_core_components as dcc
import dash_html_components as html
import dash_table as dt

from openomics_web.utils.str_utils import longest_common_prefix


def DataTableColumnSelect(columns):
    """
    Args:
        columns:
    """
    longest_common_prefixes = longest_common_prefix(columns)

    return html.Div([
        html.Div(["Select the gene id/name column to index by:"]),
        dcc.Dropdown(
            id="data-table-genes-col-name",
            options=[{
                "label": col,
                "value": col
            } for col in columns],
            style={
                "width": "100%",
            },
            value=columns[0],
        ),
        html.Div(["Select the column prefixes to import:"]),
        dcc.Dropdown(
            id="data-table-columns-select",
            options=[{
                "label": col,
                "value": col
            } for col in longest_common_prefixes],
            style={
                "width": "100%",
            },
            multi=True,
        ),
    ])


def ExpressionDataTable(df):
    """
    Args:
        df:
    """
    return html.Div(
        className="row",
        children=[
            html.Div(
                dt.DataTable(
                    id="expression-datatable",
                    columns=[{
                        "name": i,
                        "id": i
                    } for i in df.columns],
                    page_current=0,
                    page_size=20,
                    page_action="custom",
                    filter_action="custom",
                    filter_query="",
                    sort_action="custom",
                    sort_mode="multi",
                    sort_by=[],
                    style_as_list_view=True,
                    style_cell={
                        "overflow": "hidden",
                        "textOverflow": "clip",
                        "whiteSpace": "normal",
                    },
                    style_data={"width": "30px"},
                    style_data_conditional=[
                        {
                            "if": {
                                "row_index": "odd"
                            },
                            "backgroundColor": "rgb(248, 248, 248)",
                        },
                    ],
                    style_table={
                        "maxHeight": "800px",
                        "width": "800px",
                        "marginTop": "5px",
                        "marginBottom": "10px",
                        "overflowX": "scroll",
                    },
                    style_header={
                        "backgroundColor": "white",
                        "fontWeight": "bold"
                    },
                    row_selectable="multi",
                    selected_rows=[],
                    # virtualization=True,
                ),
                style={
                    "height": 750,
                    "overflowY": "scroll"
                },
                className="six columns",
            ),
            html.Div(id="table-paging-with-graph-container",
                     className="five columns"),
        ],
    )


operators = [
    ["ge ", ">="],
    ["le ", "<="],
    ["lt ", "<"],
    ["gt ", ">"],
    ["ne ", "!="],
    ["eq ", "="],
    ["contains "],
    ["datestartswith "],
]


def split_filter_part(filter_part):
    """
    Args:
        filter_part:
    """
    for operator_type in operators:
        for operator in operator_type:
            if operator in filter_part:
                name_part, value_part = filter_part.split(operator, 1)
                name = name_part[name_part.find("{") + 1:name_part.rfind("}")]

                value_part = value_part.strip()
                v0 = value_part[0]
                if v0 == value_part[-1] and v0 in ("'", '"', "`"):
                    value = value_part[1:-1].replace("\\" + v0, v0)
                else:
                    try:
                        value = float(value_part)
                    except ValueError:
                        value = value_part

                # word operators need spaces after them in the filter string,
                # but we don't want these later
                return name, operator_type[0].strip(), value

    return [None] * 3


def expression_data_view():
    return html.Div(
        id="table-container",
        children=[
            dt.DataTable(
                id="data-table",
                row_selectable="multi",
                # sorting=True,
                # filtering=True,
                css=[{
                    "selector":
                    ".dash-cell div.dash-cell-value",
                    "rule":
                    "display: inline; "
                    "white-space: inherit; "
                    "overflow: auto; "
                    "text-overflow: inherit;",
                }],
                style_cell={
                    "whiteSpace": "no-wrap",
                    "overflow": "hidden",
                    "textOverflow": "ellipsis",
                    "maxWidth": 100,
                    "fontWeight": 100,
                    "fontSize": "11pt",
                    "fontFamily": "Courier New",
                    "backgroundColor": "#1F2132",
                },
                style_header={
                    "backgroundColor": "#1F2132",
                    "textAlign": "center"
                },
                style_table={
                    "maxHeight": "310px",
                    "width": "320px",
                    "marginTop": "5px",
                    "marginBottom": "10px",
                },
                # n_fixed_rows=1,
                # n_fixed_columns=1
            )
        ],
    )
