import validators
import io
import os
from typing import List, Union

import dask.dataframe as dd
import pandas as pd
from .utils import get_pkg_data_filename

BCR_PATIENT_BARCODE_COL = "bcr_patient_barcode"
HISTOLOGIC_SUBTYPE_COL = "histologic_subtype"
PATHOLOGIC_STAGE_COL = "pathologic_stage"
TUMOR_NORMAL_COL = "tumor_normal"
PREDICTED_SUBTYPE_COL = "predicted_subtype"

TUMOR = "Tumor"
NORMAL = "Normal"


class ClinicalData:
    """This class manages the clinical data tables to handle the patient's
    phenotype data, as well as the treatment, and sample data associated to each
    patient.
    """

    pathologic_stage_map = {
        "Stage IA": "Stage I",
        "Stage IB": "Stage I",
        "Stage IIA": "Stage II",
        "Stage IIB": "Stage II",
        "Stage IIIA": "Stage III",
        "Stage IIIB": "Stage III",
    }

    def __init__(
        self,
        file_path: Union[str, io.StringIO, pd.DataFrame, dd.DataFrame],
        patient_index: str,
        columns: List[str] = None,
    ):
        """
        Args:
            file_path (str, io.StringIO, pd.DataFrame): either a path to the
                patients clinical data file, or a DataFrame.
            patient_index (str): the patient's ID column name
            columns (List[str]): default None. Specifies the columns to import,
                if None, then import all columns.
        """
        # self.cohort_name = cohort_name
        self.patient_column = patient_index

        if columns and patient_index not in columns:
            columns.append(patient_index)

        if isinstance(file_path, io.StringIO):
            file_path.seek(
                0
            )  # Needed since the file was previous read to extract columns information
            self.patient = pd.read_table(
                file_path,
                skiprows=[1, 2],
                na_values=[
                    "[Not Available]",
                    "[Unknown]",
                    "[Not Applicable]",
                    "[Discrepancy]",
                ],
                usecols=columns,
            )
        elif isinstance(file_path, str) and os.path.exists(file_path):
            self.patient = pd.read_table(
                file_path,
                skiprows=[1, 2],
                na_values=[
                    "[Not Available]",
                    "[Unknown]",
                    "[Not Applicable]",
                    "[Discrepancy]",
                ],
                usecols=columns,
            )
        elif isinstance(file_path, (pd.DataFrame, dd.DataFrame)):
            self.patient = file_path

        elif isinstance(file_path, str) and validators.url(file_path):
            dataurl, filename = os.path.split(file_path)
            file = get_pkg_data_filename(dataurl + "/", filename)
            self.patient = pd.read_table(file)

        else:
            raise IOError(file_path)

        self.patient_barcodes = self.patient[patient_index].tolist()
        self.patient.set_index(patient_index, inplace=True)

        # Rename columns
        self.patient.rename(
            {
                "ajcc_pathologic_tumor_stage": PATHOLOGIC_STAGE_COL,
                "histological_type": HISTOLOGIC_SUBTYPE_COL,
                "histologic_diagnosis.1": HISTOLOGIC_SUBTYPE_COL,
            },
            axis=1,
            inplace=True,
        )

        self.patient.replace(
            {PATHOLOGIC_STAGE_COL: ClinicalData.pathologic_stage_map},
            inplace=True)

    @classmethod
    def name(self):
        """Returns the name of the class, i.e. 'ClinicalData'"""
        return self.__class__.__name__

    def build_clinical_samples(self, all_samples, index="bcr_patient_barcode"):
        """Build table with samples clinical data from patients :param
        all_samples:

        Args:
            all_samples:
            index:
        """
        self.samples = pd.DataFrame(index=all_samples)
        self.samples.index.name = index
        self.samples.index = self.samples.index.str[:
                                                    -4]  # Cut sample barcode for TCGA

        num_samples = self.samples.shape[0]

        # Merge patients clinical data with patient barcode as index
        # target = pd.merge(target, self.patient,
        #                      how="left", left_on="patient_barcode", right_on="patient_barcode")

        self.samples = self.samples.join(self.patient,
                                         on=index,
                                         how="left",
                                         rsuffix="_")
        if self.samples.shape[0] != num_samples:
            raise Exception(
                "Clinical data merging has wrong number of samples")

        # self.samples.dropna(axis=0, subset=["bcr_patient_barcode"], inplace=True) # Remove samples without clinical data

        self.samples = self.samples[
            self.samples[PATHOLOGIC_STAGE_COL] != "[Discrepancy]"]
        self.samples.loc[
            self.samples.index.str.contains("-11"),
            TUMOR_NORMAL_COL] = NORMAL  # Change stage label of normal samples to "Normal"
        self.samples.loc[
            self.samples.index.str.contains("-01"),
            TUMOR_NORMAL_COL] = TUMOR  # Change stage label of normal samples to "Normal"

    def add_drug_response_data(
        self,
        file_path="nationwidechildrens.org_clinical_drug.txt",
        patient_column="bcr_patient_barcode",
        columns=None,
        drug_name_col=None,
        response_column=None,
    ):
        """
        Args:
            file_path:
            patient_column:
            columns:
            drug_name_col:
            response_column:
        """
        if columns is None:
            columns = [
                "bcr_patient_barcode",
                "pharmaceutical_therapy_drug_name",
                "pharmaceutical_therapy_type",
                "treatment_best_response",
            ]

        if not os.path.exists(file_path):
            raise FileNotFoundError(file_path)

        self.drug_name_col = drug_name_col
        self.response_column = response_column

        self.drugs = pd.read_table(
            file_path,
            sep="\t",
            skiprows=[1, 2],
            na_values=["[Not Available]", "[Unknown]", "[Not Applicable]"],
            usecols=columns,
        )
        self.drugs.set_index(patient_column, inplace=True)

    def add_biospecimen_data(
        self,
        file_path="genome.wustl.edu_biospecimen_sample.txt",
        patient_col_name="bcr_patient_barcode",
        columns=["bcr_sample_barcode", "sample_type"],
    ):
        """
        Args:
            file_path:
            patient_col_name:
            columns:
        """
        if not os.path.exists(file_path):
            raise FileNotFoundError(file_path)

        self.biospecimen = pd.read_table(
            file_path,
            sep="\t",
            skiprows=[
                1,
            ],
            na_values=["[Not Available]", "[Unknown]", "[Not Applicable]"],
            usecols=columns,
        )
        self.sample_barcodes = self.biospecimen[patient_col_name].tolist()
        self.biospecimen.set_index(patient_col_name, inplace=True)

    def get_patient_barcodes(self):
        return self.patient_barcodes

    def get_sample_barcodes(self):
        return self.sample_barcodes


# class DrugResponse():
#     def __init__(self, drugs_file_path="nationwidechildrens.org_clinical_drug.txt", patient_column="bcr_patient_barcode",
#                  columns=['bcr_patient_barcode', 'pharmaceutical_therapy_drug_name', 'pharmaceutical_therapy_type', 'treatment_best_response'],
#                  drug_name_col=None, response_column=None):
#         self.drug_name_col = drug_name_col
#         self.response_column = response_column
#
#         self.drugs = pd.read_table(drugs_file_path,
#                                    sep="\t",
#                                    skiprows=[1, 2],
#                                    na_values=["[Not Available]", "[Unknown]", "[Not Applicable]"],
#                                    usecols=columns
#                                    )
#         self.drugs.set_index(patient_column, inplace=True)

# class Biospecimen():
#     def __init__(self, biospecimens_file="genome.wustl.edu_biospecimen_sample.txt", patient_col_name="bcr_patient_barcode",
#                  columns=['bcr_sample_barcode', 'sample_type']):
#         self.biospecimen = pd.read_table(biospecimens_file, sep="\t", skiprows=[1, ],
#                                          na_values=["[Not Available]", "[Unknown]", "[Not Applicable]"],
#                                          usecols=columns
#                                          )
#         self.sample_barcodes = self.biospecimen[patient_col_name].tolist()
#         self.biospecimen.set_index(patient_col_name, inplace=True)
