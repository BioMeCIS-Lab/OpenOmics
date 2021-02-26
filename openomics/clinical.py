import io
import os

import dask.dataframe as dd
from openomics import backend as pd

TUMOR = "Tumor"
NORMAL = "Normal"

BCR_PATIENT_BARCODE = "bcr_patient_barcode"
HISTOLOGIC_SUBTYPE = "histologic_subtype"
PATHOLOGIC_STAGE = "pathologic_stage"
TUMOR_NORMAL = "tumor_normal"
PREDICTED_SUBTYPE = "predicted_subtype"


class ClinicalData:
    pathologic_stage_map = {
        "Stage IA": "Stage I",
        "Stage IB": "Stage I",
        "Stage IIA": "Stage II",
        "Stage IIB": "Stage II",
        "Stage IIIA": "Stage III",
        "Stage IIIB": "Stage III",
    }

    def __init__(self, file_path, patient_index, columns=None):
        """This class manages the clinical data tables to handle the phenotype,
        treatment, and sample data associated to a patient.

        Args:
            file_path (str, io.StringIO, pd.DataFrame): either a path to the patients clinical data file, or a DataFrame.
            patient_index (str): the patient's ID column name
            columns (list): default None. Specifies the columns to import, if
                None, then import all columns.
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
        elif type(file_path) == str and os.path.exists(file_path):
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

        else:
            raise IOError(file_path)

        self.patient_barcodes = self.patient[patient_index].tolist()
        self.patient.set_index(patient_index, inplace=True)

        # Rename columns
        self.patient.rename(
            {
                "ajcc_pathologic_tumor_stage": PATHOLOGIC_STAGE,
                "histological_type": HISTOLOGIC_SUBTYPE,
                "histologic_diagnosis.1": HISTOLOGIC_SUBTYPE,
            },
            axis=1,
            inplace=True,
        )

        self.patient.replace(
            {PATHOLOGIC_STAGE: ClinicalData.pathologic_stage_map},
            inplace=True)

    @classmethod
    def name(self):
        return self.__class__.__name__

    def build_clinical_samples(self, all_samples, index="bcr_patient_barcode"):
        # Build table with samples clinical data from patients
        """
        Args:
            all_samples:
            index:
        """
        self.samples = pd.DataFrame(index=all_samples)
        self.samples.index.name = index
        self.samples.index = self.samples.index.str[:
                                                    -4]  # Cut sample barcode for TCGA

        no_samples = self.samples.shape[0]

        # Merge patients clinical data with patient barcode as index
        # target = pd.merge(target, self.patient,
        #                      how="left", left_on="patient_barcode", right_on="patient_barcode")
        self.samples = self.samples.join(self.patient,
                                         on=index,
                                         how="left",
                                         rsuffix="_")
        if self.samples.shape[0] != no_samples:
            raise Exception(
                "Clinical data merging has wrong number of samples")

        # self.samples.dropna(axis=0, subset=["bcr_patient_barcode"], inplace=True) # Remove samples without clinical data

        self.samples = self.samples[
            self.samples[PATHOLOGIC_STAGE] != "[Discrepancy]"]
        self.samples.loc[
            self.samples.index.str.contains("-11"),
            TUMOR_NORMAL] = NORMAL  # Change stage label of normal samples to "Normal"
        self.samples.loc[
            self.samples.index.str.contains("-01"),
            TUMOR_NORMAL] = TUMOR  # Change stage label of normal samples to "Normal"

    def add_drug_response_data(
        self,
        file_path="nationwidechildrens.org_clinical_drug.txt",
        patient_column="bcr_patient_barcode",
        columns=[
            "bcr_patient_barcode",
            "pharmaceutical_therapy_drug_name",
            "pharmaceutical_therapy_type",
            "treatment_best_response",
        ],
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
