import dash
import dash_html_components as html
from dash.dependencies import Input, Output, State

from openomics import MultiOmics
from openomics_web.layouts import app_layout
from openomics_web.layouts.clinical_view import (
    ClinicalDataColumnSelect,
    ClinicalDataTable,
)
from openomics_web.layouts.datatable_view import (
    ExpressionDataTable,
    DataTableColumnSelect,
    split_filter_part,
)
from openomics_web.server import server
from openomics_web.utils.io import (
    get_table_columns,
    get_expression_data,
    get_clinical_data,
)

external_stylesheets = ["https://codepen.io/chriddyp/pen/bWLwgP.css"]

# running directly with Python
app = dash.Dash(__name__,
                server=server,
                external_stylesheets=external_stylesheets)

app.layout = app_layout.app_main()

user_multiomics = MultiOmics(cohort_name="TEST", )


@app.callback(
    [
        Output("data-table-column-select", "children"),
        Output("upload-data-table", "children"),
    ],
    [
        Input("upload-data-table", "contents"),
        Input("upload-data-table", "filename")
    ],
    [
        State("data-table-type", "value"),
    ],
)
def update_datatable_metadata(
    list_of_contents,
    list_of_names,
    data_type,
):
    """
    Args:
        list_of_contents:
        list_of_names:
        data_type:
    """
    if list_of_contents is None:
        return None, ["Drag and Drop or ", html.A("Select Files")]

    try:
        columns = get_table_columns(list_of_contents, list_of_names)

    except Exception as e:
        print(e)
        return None, "There was an error processing this file."

    return DataTableColumnSelect(columns), "Uploaded {}".format(list_of_names)


@app.callback(
    Output("output-data-upload", "children"),
    [Input("upload-data-table-submit", "n_clicks")],
    [
        State("data-table-cohort", "value"),
        State("data-table-type", "value"),
        State("upload-data-table", "contents"),
        State("upload-data-table", "filename"),
        State("data-table-genes-col-name", "value"),
        State("data-table-columns-select", "value"),
        State("data-table-transpose", "value"),
    ],
)
def import_datatable_upload(
    n_clicks,
    cohort_name,
    data_type,
    list_of_contents,
    list_of_names,
    genes_col_name,
    columns_select,
    transposed,
):
    """
    Args:
        n_clicks:
        cohort_name:
        data_type:
        list_of_contents:
        list_of_names:
        genes_col_name:
        columns_select:
        transposed:
    """
    if list_of_contents is None:
        return []
    try:

        omics_data = get_expression_data(
            list_of_contents,
            list_of_names,
            data_type,
            cohort_name,
            genes_col_name,
            columns_select,
            transposed,
        )
        user_multiomics.add_omic(omics_data)
    except Exception as e:
        print(e)
        return html.Div(["There was an error processing this file."])

    return ExpressionDataTable(omics_data.expressions.head(20))


@app.callback(
    Output("expression-datatable", "data"),
    [
        Input("expression-datatable", "page_current"),
        Input("expression-datatable", "page_size"),
        Input("expression-datatable", "sort_by"),
        Input("expression-datatable", "filter_query"),
    ],
)
def update_table(page_current, page_size, sort_by, filter):
    """
    Args:
        page_current:
        page_size:
        sort_by:
        filter:
    """
    filtering_expressions = filter.split(" && ")
    print(user_multiomics.get_omics_list())
    dff = user_multiomics[user_multiomics.get_omics_list()[0]]
    for filter_part in filtering_expressions:
        col_name, operator, filter_value = split_filter_part(filter_part)

        if operator in ("eq", "ne", "lt", "le", "gt", "ge"):
            # these operators match pandas series operator method names
            dff = dff.loc[getattr(dff[col_name], operator)(filter_value)]
        elif operator == "contains":
            dff = dff.loc[dff[col_name].str.contains(filter_value)]
        elif operator == "datestartswith":
            # this is a simplification of the front-end filtering logic,
            # only works with complete fields in standard format
            dff = dff.loc[dff[col_name].str.startswith(filter_value)]

    if sort_by:
        dff = dff.sort_values(
            [col["column_id"] for col in sort_by],
            ascending=[col["direction"] == "asc" for col in sort_by],
            inplace=False,
        )

    return dff.iloc[page_current * page_size:(page_current + 1) *
                    page_size].to_dict("records")


@app.callback(
    [
        Output("clinical-column-select", "children"),
        Output("upload-clinical", "children"),
    ],
    [
        Input("upload-clinical", "contents"),
        Input("upload-clinical", "filename")
    ],
)
def update_clinical_upload_metadata(
    file_content,
    file_name,
):
    """
    Args:
        file_content:
        file_name:
    """
    if file_content is None:
        return None, ["Drag and Drop or ", html.A("Select Files")]

    try:
        columns = get_table_columns(
            [
                file_content,
            ],
            [
                file_name,
            ],
        )

    except Exception as e:
        print(e)
        return None, "There was an error processing this file."

    return ClinicalDataColumnSelect(columns), "Uploaded {}".format(file_name)


@app.callback(
    Output("output-clinical-upload", "children"),
    [Input("clinical-submit-button", "n_clicks")],
    [
        State("clinical-cohort", "value"),
        State("clinical-data-type", "value"),
        State("upload-clinical", "contents"),
        State("upload-clinical", "filename"),
        State("clinical-patient-col-name", "value"),
        State("clinical-data-columns-select", "value"),
    ],
)
def import_datatable_upload(
    n_clicks,
    cohort_name,
    data_type,
    list_of_contents,
    list_of_names,
    patient_id_col,
    columns_select,
):
    """
    Args:
        n_clicks:
        cohort_name:
        data_type:
        list_of_contents:
        list_of_names:
        patient_id_col:
        columns_select:
    """
    if list_of_contents is None:
        return []
    try:
        clinical_data = get_clinical_data(
            list_of_contents,
            list_of_names,
            data_type,
            cohort_name,
            patient_id_col,
            columns_select,
        )
        user_multiomics.add_clinical_data(clinical_data)
    except Exception as e:
        print(e)
        return html.Div(["There was an error processing this file."])

    return ClinicalDataTable(clinical_data.patient.head(20))


if __name__ == "__main__":
    app.run_server(debug=True, port=8050)
