import os

import pytest

from openomics import ClinicalData

cohort_folder_path = "tests/data/TCGA_LUAD"


@pytest.fixture
def generate_TCGA_clinical():
    clinical = ClinicalData(
        file_path=os.path.join(
            cohort_folder_path,
            "nationwidechildrens.org_clinical_patient_luad.txt"),
        patient_index="bcr_patient_barcode",
    )
    return clinical
