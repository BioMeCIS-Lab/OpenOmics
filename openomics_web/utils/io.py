import base64
import io
import os

import pandas as pd

from openomics import ClinicalData
from openomics import MicroRNA, MessengerRNA, LncRNA, Protein


def get_table_columns(list_of_contents, list_of_names):
    """
    Args:
        list_of_contents:
        list_of_names:
    """
    content = list_of_contents[0]
    filename = list_of_names[0]
    _, file_extension = os.path.splitext(filename)
    content_type, content_string = content.split(",")
    decoded = base64.b64decode(content_string)

    if file_extension == ".tsv":
        first_line = decoded.decode("utf-8").partition("\n")[0]
        columns = first_line.split("\t")
    elif file_extension == ".csv":
        first_line = decoded.decode("utf-8").partition("\n")[0]
        columns = first_line.split(",")
    elif file_extension == ".xls" or file_extension == ".xlsx":
        columns = pd.read_excel(io.BytesIO(decoded),
                                low_memory=True).columns.tolist()
    else:
        columns = pd.read_table(io.StringIO(decoded.decode("utf-8")),
                                low_memory=True).columns.tolist()
    return columns


def get_expression_data(
    list_of_contents,
    list_of_names,
    data_type,
    cohort_name=None,
    genes_col_name=None,
    columns_selected=None,
    transposed=None,
):
    """
    Args:
        list_of_contents:
        list_of_names:
        data_type:
        cohort_name:
        genes_col_name:
        columns_selected:
        transposed:
    """
    if columns_selected:
        columns = "|".join(columns_selected)
    else:
        columns = None

    if transposed == "True":
        transposed = True
    elif transposed == "False":
        transposed = False

    file = handle_filestreams(list_of_contents, list_of_names)

    if data_type == MicroRNA.name():
        expression_data = MicroRNA(file, transpose=transposed, usecols=columns)
    elif data_type == MessengerRNA.name():
        expression_data = MessengerRNA(file,
                                       transpose=transposed,
                                       usecols=columns)
    elif data_type == LncRNA.name():
        expression_data = LncRNA(file, transpose=transposed, usecols=columns)
    elif data_type == Protein.name():
        expression_data = Protein(file, transpose=transposed, usecols=columns)
    else:
        print(data_type)

    return expression_data


def get_clinical_data(
    file_content,
    file_name,
    data_type,
    cohort_name,
    patient_id_col=None,
    columns_selected=None,
):
    """
    Args:
        file_content:
        file_name:
        data_type:
        cohort_name:
        patient_id_col:
        columns_selected:
    """
    file = handle_filestreams(
        [
            file_content,
        ],
        [
            file_name,
        ],
    )
    clinical_data = ClinicalData(file,
                                 patient_index=patient_id_col,
                                 columns=columns_selected)

    return clinical_data


def handle_filestreams(list_of_contents, list_of_names):
    """
    Args:
        list_of_contents:
        list_of_names:
    """
    if len(list_of_contents) == 1:
        content = list_of_contents[0]
        filename = list_of_names[0]
    else:
        raise Exception("Multiple files not supported")  # TODO

    content_type, content_string = content.split(",")
    decoded = base64.b64decode(content_string)
    if "csv" in filename:  # Assume that the user uploaded a CSV file
        file = io.StringIO(decoded.decode("utf-8"))
    elif "xls" in filename:  # Assume that the user uploaded an excel file
        file = io.BytesIO(decoded)
    elif "tsv" in filename:  # Assume that the user uploaded an tsv file
        file = io.StringIO(decoded.decode("utf-8"))
    elif "txt" in filename:  # Assume that the user uploaded either a tsv or csv file
        file = io.StringIO(decoded.decode("utf-8"))
    else:
        raise IOError("Unable to read table file.")
    return file
