from openomics.database import DisGeNet, HMDD, LncRNADisease, MalaCards
from .test_multiomics import *


@pytest.fixture
def generate_DisGeNet_ftp():
    return DisGeNet(
        path="https://www.disgenet.org/static/disgenet_ap1/files/downloads/",
        curated=True,
    )


@pytest.fixture
def generate_LncRNADisease_ftp():
    return LncRNADisease(path="http://www.cuilab.cn/files/images/ldd/",
                         species="Human")


@pytest.fixture
def generate_HMDD_ftp():
    return HMDD(path="http://www.cuilab.cn/static/hmdd3/data/")


@pytest.fixture
def generate_MalaCards_ftp():
    return MalaCards()


def test_import_DisGeNet(generate_DisGeNet_ftp):
    """
    Args:
        generate_DisGeNet_ftp:
    """
    assert (generate_DisGeNet_ftp.data_path ==
            "https://www.disgenet.org/static/disgenet_ap1/files/downloads/")
    assert not generate_DisGeNet_ftp.get_disease_assocs(
        index="gene_name").empty


def test_annotate_DisGeNet(generate_TCGA_LUAD, generate_DisGeNet_ftp):
    """
    Args:
        generate_TCGA_LUAD:
        generate_DisGeNet_ftp:
    """
    generate_TCGA_LUAD.MessengerRNA.annotate_diseases(
        generate_DisGeNet_ftp,
        index="gene_name",
    )
    assert {"disease_associations"
            }.issubset(generate_TCGA_LUAD.MessengerRNA.annotations.columns)


def test_import_HMDD(generate_HMDD_ftp):
    """
    Args:
        generate_HMDD_ftp:
    """
    assert generate_HMDD_ftp.data_path == "http://www.cuilab.cn/static/hmdd3/data/"
    assert not generate_HMDD_ftp.get_disease_assocs(index="gene_name").empty


def test_annotate_HMDD(generate_TCGA_LUAD, generate_HMDD_ftp):
    """
    Args:
        generate_TCGA_LUAD:
        generate_HMDD_ftp:
    """
    generate_TCGA_LUAD.MicroRNA.annotate_diseases(
        generate_HMDD_ftp,
        index="gene_name",
    )
    assert {"disease_associations"
            }.issubset(generate_TCGA_LUAD.MicroRNA.annotations.columns)


def test_LncRNADisease(generate_LncRNADisease_ftp):
    """
    Args:
        generate_LncRNADisease_ftp:
    """
    assert (generate_LncRNADisease_ftp.data_path ==
            "http://www.cuilab.cn/files/images/ldd/")
    assert not generate_LncRNADisease_ftp.get_disease_assocs(
        index="gene_name").empty


def test_annotate_LncRNADisease(generate_TCGA_LUAD,
                                generate_LncRNADisease_ftp):
    """
    Args:
        generate_TCGA_LUAD:
        generate_LncRNADisease_ftp:
    """
    generate_TCGA_LUAD.LncRNA.annotate_diseases(
        generate_LncRNADisease_ftp,
        index="gene_name",
    )
    assert {"disease_associations"
            }.issubset(generate_TCGA_LUAD.LncRNA.annotations.columns)


def test_import_MalaCards(generate_MalaCards_ftp):
    """
    Args:
        generate_MalaCards_ftp:
    """
    assert generate_MalaCards_ftp.data_path == "http://zdzlab.einstein.yu.edu/1/hedd/"
