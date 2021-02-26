from openomics.database import GENCODE, MirBase
from .test_multiomics import *


@pytest.fixture
def generate_GENCODE():
    return GENCODE(
        path=
        "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/",
        file_resources={
            "long_noncoding_RNAs.gtf":
            "gencode.v32.long_noncoding_RNAs.gtf.gz",
            "lncRNA_transcripts.fa": "gencode.v32.lncRNA_transcripts.fa.gz",
            "transcripts.fa": "gencode.v32.transcripts.fa.gz",
        },
    )


@pytest.fixture
def generate_GENCODE_dask():
    return GENCODE(
        path=
        "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/",
        file_resources={
            "long_noncoding_RNAs.gtf":
            "gencode.v32.long_noncoding_RNAs.gtf.gz",
            "lncRNA_transcripts.fa": "gencode.v32.lncRNA_transcripts.fa.gz",
            "transcripts.fa": "gencode.v32.transcripts.fa.gz",
        },
        npartitions=8,
    )


def test_import_GENCODE(generate_GENCODE):
    """
    Args:
        generate_GENCODE:
    """
    assert (
        generate_GENCODE.data_path ==
        "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/")


def test_import_GENCODE_dask(generate_GENCODE_dask):
    assert (
        generate_GENCODE_dask.data_path ==
        "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/")


def test_annotate_GENCODE(generate_TCGA_LUAD, generate_GENCODE):
    """
    Args:
        generate_TCGA_LUAD:
        generate_GENCODE:
    """
    generate_TCGA_LUAD.LncRNA.annotate_genomics(
        generate_GENCODE,
        index="gene_id",
        columns=["feature", "start", "end", "strand", "tag", "havana_gene"],
    )
    assert {"feature", "start", "end", "strand", "tag", "havana_gene"
            }.issubset(generate_TCGA_LUAD.LncRNA.get_annotations().columns)


@pytest.fixture
def generate_MirBase_ftp():
    return MirBase(path="ftp://mirbase.org/pub/mirbase/CURRENT/")


def test_import_mirbase_db(generate_MirBase_ftp):
    """
    Args:
        generate_MirBase_ftp:
    """
    assert generate_MirBase_ftp.data_path == "ftp://mirbase.org/pub/mirbase/CURRENT/"
