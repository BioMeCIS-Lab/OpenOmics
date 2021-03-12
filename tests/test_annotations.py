from openomics.database import RNAcentral, GTEx, GeneOntology
from .test_multiomics import *


@pytest.fixture
def generate_RNACentral_ftp():
    rnacentral = RNAcentral(
        path="ftp://ftp.ebi.ac.uk/pub/databases/RNAcentral/current_release/")
    rnacentral.data = rnacentral.data.sample(frac=0.01)
    return rnacentral


@pytest.fixture
def generate_GTEx_expressions():
    return GTEx(
        path="https://storage.googleapis.com/gtex_analysis_v8/rna_seq_data/", )


@pytest.fixture
def generate_GeneOntology():
    go = GeneOntology(path="http://geneontology.org/gene-associations/")
    go.data = go.data.sample(frac=0.01)
    return go


def test_import_rnacentral_db(generate_RNACentral_ftp):
    """
    Args:
        generate_RNACentral_ftp:
    """
    assert (generate_RNACentral_ftp.data_path ==
            "ftp://ftp.ebi.ac.uk/pub/databases/RNAcentral/current_release/")


def test_import_GTEx(generate_GTEx_expressions):
    """
    Args:
        generate_GTEx_expressions:
    """
    assert (generate_GTEx_expressions.data_path ==
            "https://storage.googleapis.com/gtex_analysis_v8/rna_seq_data/")
    assert not generate_GTEx_expressions.get_expressions(
        index="gene_name").empty
    assert not generate_GTEx_expressions.get_expressions(index="gene_id").empty


def test_import_GeneOntology(generate_GeneOntology):
    """
    Args:
        generate_GeneOntology:
    """
    assert (generate_GeneOntology.data_path ==
            "http://geneontology.org/gene-associations/")


def test_annotate_rnacentral(generate_TCGA_LUAD, generate_RNACentral_ftp):
    """
    Args:
        generate_TCGA_LUAD:
        generate_RNACentral_ftp:
    """
    generate_TCGA_LUAD.MicroRNA.annotate_attributes(
        database=generate_RNACentral_ftp,
        on="gene_name",
        columns=["gene_name", "RNA type"],
    )
    assert {"RNA type"
            }.issubset(generate_TCGA_LUAD.MicroRNA.annotations.columns)


def test_annotate_expressions_GTEx(generate_TCGA_LUAD,
                                   generate_GTEx_expressions):
    """
    Args:
        generate_TCGA_LUAD:
        generate_GTEx_expressions:
    """
    generate_TCGA_LUAD.LncRNA.annotate_expressions(
        database=generate_GTEx_expressions, index="gene_id")
    assert not generate_TCGA_LUAD.LncRNA.get_annotation_expressions().empty


def test_annotate_GeneOntology(generate_TCGA_LUAD, generate_GeneOntology):
    """
    Args:
        generate_TCGA_LUAD:
        generate_GeneOntology:
    """
    generate_TCGA_LUAD.MessengerRNA.annotate_attributes(
        database=generate_GeneOntology, on="gene_name", columns=["go_id"])
    assert {"go_id"
            }.issubset(generate_TCGA_LUAD.MessengerRNA.annotations.columns)
    assert not generate_TCGA_LUAD.MessengerRNA.annotations["go_id"].empty
