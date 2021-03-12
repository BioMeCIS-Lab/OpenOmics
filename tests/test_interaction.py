from openomics.database import MiRTarBase, STRING, LncRNA2Target
from .test_multiomics import *


@pytest.fixture
def generate_LncRNA2Target():
    return LncRNA2Target(path="http://123.59.132.21/lncrna2target/data/",
                         version="low_throughput")


def test_import_LncRNA2Target(generate_LncRNA2Target):
    """
    Args:
        generate_LncRNA2Target:
    """
    assert generate_LncRNA2Target.data_path is not None


@pytest.fixture
def generate_STRING():
    string = STRING(edge_attr=["score"])
    string.data = string.data.sample(frac=0.01)
    return string


@pytest.fixture
def generate_MiRTarBase():
    return MiRTarBase(
        path=
        "/data/datasets/Bioinformatics_ExternalData/miRTarBase/",  # Hard-coded
        strip_mirna_name=True,
        filters={"Species (Target Gene)": "Homo sapiens"},
    )


# Test disabled since obtaining MiRTarBase via ftp is unreachable
# def test_import_MiRTarBase(generate_MiRTarBase):
#     """
#     Args:
#         generate_MiRTarBase:
#     """
#     assert generate_MiRTarBase.data_path is not None


def test_import_STRING(generate_STRING):
    """
    Args:
        generate_STRING:
    """
    assert generate_STRING.data_path is not None


def test_annotate_STRING(generate_TCGA_LUAD, generate_STRING):
    """
    Args:
        generate_TCGA_LUAD:
        generate_STRING:
    """
    generate_TCGA_LUAD.Protein.annotate_sequences(generate_STRING,
                                                  index="protein_name")
    assert not generate_TCGA_LUAD.Protein.annotations["sequence"].empty


def test_get_interactions_lnc2target(generate_LncRNA2Target):
    assert generate_LncRNA2Target.get_interactions() is not None
