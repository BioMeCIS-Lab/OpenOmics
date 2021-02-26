from openomics.database import RNAcentral, GTEx, GeneOntology
from .test_multiomics import *


@pytest.fixture
def generate_RNACentral_ftp():
    return RNAcentral(
        path="ftp://ftp.ebi.ac.uk/pub/databases/RNAcentral/current_release/")


def test_import_rnacentral_db(generate_RNACentral_ftp):
    """
    Args:
        generate_RNACentral_ftp:
    """
    assert (generate_RNACentral_ftp.data_path ==
            "ftp://ftp.ebi.ac.uk/pub/databases/RNAcentral/current_release/")


def test_rnacentral_annotate(generate_TCGA_LUAD, generate_RNACentral_ftp):
    """
    Args:
        generate_TCGA_LUAD:
        generate_RNACentral_ftp:
    """
    generate_TCGA_LUAD.MessengerRNA.annotate_genomics(
        database=generate_RNACentral_ftp,
        index="gene_name",
        columns=["gene_name", "transcript_id", "RNA type", "go_id", "Rfams"],
    )
    assert {"transcript_id", "RNA type", "go_id", "Rfams"}.issubset(
        generate_TCGA_LUAD.MessengerRNA.get_annotations().columns)


@pytest.fixture
def generate_GTEx_expressions():
    return GTEx(
        path="https://storage.googleapis.com/gtex_analysis_v8/rna_seq_data/", )


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


def test_GTEx_annotate(generate_TCGA_LUAD, generate_GTEx_expressions):
    """
    Args:
        generate_TCGA_LUAD:
        generate_GTEx_expressions:
    """
    generate_TCGA_LUAD.MessengerRNA.annotate_expressions(
        database=generate_GTEx_expressions, index="gene_name")
    generate_TCGA_LUAD.LncRNA.annotate_expressions(
        database=generate_GTEx_expressions, index="gene_id")
    generate_TCGA_LUAD.MicroRNA.annotate_expressions(
        database=generate_GTEx_expressions, index="gene_name")
    assert not generate_TCGA_LUAD.MessengerRNA.get_annotation_expressions(
    ).empty
    assert not generate_TCGA_LUAD.LncRNA.get_annotation_expressions().empty
    assert not generate_TCGA_LUAD.MicroRNA.get_annotation_expressions().empty


@pytest.fixture
def generate_GeneOntology():
    return GeneOntology(path="http://geneontology.org/gene-associations/")


def test_import_GeneOntology(generate_GeneOntology):
    """
    Args:
        generate_GeneOntology:
    """
    assert (generate_GeneOntology.data_path ==
            "http://geneontology.org/gene-associations/")


def test_annotate_GeneOntology(generate_TCGA_LUAD, generate_GeneOntology):
    """
    Args:
        generate_TCGA_LUAD:
        generate_GeneOntology:
    """
    generate_TCGA_LUAD.MessengerRNA.annotate_genomics(
        database=generate_GeneOntology, index="gene_name", columns=["go_id"])
    generate_TCGA_LUAD.LncRNA.annotate_genomics(database=generate_GeneOntology,
                                                index="gene_id",
                                                columns=["go_id"])
    generate_TCGA_LUAD.MicroRNA.annotate_genomics(
        database=generate_GeneOntology, index="gene_name", columns=["go_id"])
    assert {"go_id"}.issubset(
        generate_TCGA_LUAD.MessengerRNA.get_annotations().columns)
    assert {"go_id"
            }.issubset(generate_TCGA_LUAD.LncRNA.get_annotations().columns)
    assert {"go_id"
            }.issubset(generate_TCGA_LUAD.MicroRNA.get_annotations().columns)
    assert not generate_TCGA_LUAD.LncRNA.annotations["go_id"].empty
    assert not generate_TCGA_LUAD.MessengerRNA.annotations["go_id"].empty
    assert not generate_TCGA_LUAD.MicroRNA.annotations["go_id"].empty
