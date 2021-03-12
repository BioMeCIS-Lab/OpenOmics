"""Tests for `openomics` package."""

from openomics import MessengerRNA, MicroRNA, LncRNA, Protein, SomaticMutation
from openomics import MultiOmics

from .test_clinical import *

cohort_folder_path = "tests/data/TCGA_LUAD"


@pytest.fixture
def generate_TCGA_LUAD_MessengerRNA():
    data = MessengerRNA(
        data=os.path.join(cohort_folder_path, "LUAD__geneExp.txt"),
        transpose=True,
        usecols="GeneSymbol|TCGA",
        gene_index="GeneSymbol",
        gene_level="gene_name",
    )
    data.drop_genes(data.expressions.columns[50:])
    data.drop_samples(data.expressions.index[:100])
    return data


@pytest.fixture
def generate_TCGA_LUAD_MessengerRNA_dask():
    data = MessengerRNA(
        data=os.path.join(cohort_folder_path, "LUAD__geneExp.txt"),
        transpose=True,
        usecols="GeneSymbol|TCGA",
        gene_index="GeneSymbol",
        gene_level="gene_name",
        npartitions=4,
    )
    data.drop_genes(data.expressions.columns[50:])
    data.initialize_annotations(index="gene_name", gene_list=None)

    return data


@pytest.fixture
def generate_TCGA_LUAD_MicroRNA():
    data = MicroRNA(
        data=os.path.join(cohort_folder_path, "LUAD__miRNAExp__RPM.txt"),
        transpose=True,
        usecols="GeneSymbol|TCGA",
        gene_index="GeneSymbol",
        gene_level="gene_name",
    )
    data.drop_genes(data.expressions.columns[50:])
    data.drop_samples(data.expressions.index[:100])
    return data


@pytest.fixture
def generate_TCGA_LUAD_LncRNA():
    data = LncRNA(
        data=os.path.join(cohort_folder_path, "TCGA-rnaexpr.tsv"),
        transpose=True,
        usecols="Gene_ID|TCGA",
        gene_index="Gene_ID",
        gene_level="gene_id",
    )
    data.drop_genes(data.expressions.columns[50:])
    data.drop_samples(data.expressions.index[:100])
    return data


@pytest.fixture
def generate_TCGA_LUAD_SomaticMutation():
    data = SomaticMutation(
        data=os.path.join(cohort_folder_path,
                          "LUAD__somaticMutation_geneLevel.txt"),
        transpose=True,
        usecols="GeneSymbol|TCGA",
        gene_index="gene_name",
    )
    data.drop_genes(data.expressions.columns[50:])
    data.drop_samples(data.expressions.index[:100])
    return data


@pytest.fixture
def generate_TCGA_LUAD_Protein():
    data = Protein(
        data=os.path.join(cohort_folder_path, "protein_RPPA.txt"),
        transpose=True,
        usecols="GeneSymbol|TCGA",
        gene_index="GeneSymbol",
        gene_level="protein_name",
    )
    data.drop_genes(data.expressions.columns[50:])
    data.drop_samples(data.expressions.index[:100])
    return data


def test_import_MessengerRNA_Dask(generate_TCGA_LUAD_MessengerRNA_dask):
    """
    Args:
        generate_TCGA_LUAD_MessengerRNA_dask:
    """
    assert generate_TCGA_LUAD_MessengerRNA_dask.expressions is not None


def test_import_expression_table_size(generate_TCGA_LUAD_MessengerRNA,
                                      generate_TCGA_clinical):
    """
    Args:
        generate_TCGA_LUAD_MessengerRNA:
    """
    cohort_name = "LUAD"
    luad_data = MultiOmics(cohort_name)
    luad_data.add_clinical_data(generate_TCGA_clinical)
    luad_data.add_omic(generate_TCGA_LUAD_MessengerRNA)
    luad_data.build_samples()
    print(luad_data.data.keys())
    assert (luad_data.data[MessengerRNA.name()].shape ==
            generate_TCGA_LUAD_MessengerRNA.expressions.shape)


@pytest.fixture
def generate_TCGA_LUAD(
    generate_TCGA_clinical,
    generate_TCGA_LUAD_MessengerRNA,
    generate_TCGA_LUAD_MicroRNA,
    generate_TCGA_LUAD_LncRNA,
    generate_TCGA_LUAD_Protein,
):
    """
    Args:
        generate_TCGA_LUAD_MessengerRNA:
        generate_TCGA_LUAD_MicroRNA:
        generate_TCGA_LUAD_LncRNA:
        generate_TCGA_LUAD_Protein:
    """
    cohort_name = "LUAD"
    luad_data = MultiOmics(cohort_name)
    luad_data.add_clinical_data(generate_TCGA_clinical)
    luad_data.add_omic(generate_TCGA_LUAD_MessengerRNA)
    luad_data.add_omic(generate_TCGA_LUAD_MicroRNA)
    luad_data.add_omic(generate_TCGA_LUAD_LncRNA)
    luad_data.add_omic(generate_TCGA_LUAD_Protein)
    return luad_data


def test_TCGA_LUAD_multiomics_transcriptomics(generate_TCGA_LUAD):
    """
    Args:
        generate_TCGA_LUAD:
    """
    assert all(elem in generate_TCGA_LUAD.get_omics_list() for elem in [
        MessengerRNA.name(),
        MicroRNA.name(),
        LncRNA.name(),
        Protein.name(),
    ])
