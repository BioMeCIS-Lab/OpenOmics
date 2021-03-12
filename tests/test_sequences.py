import dask.dataframe as dd

from openomics.database import GENCODE, MirBase
from .test_multiomics import *


@pytest.fixture
def generate_GENCODE():
    gencode = GENCODE(
        path=
        "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/",
        file_resources={
            "long_noncoding_RNAs.gtf":
            "gencode.v32.long_noncoding_RNAs.gtf.gz",
            "lncRNA_transcripts.fa": "gencode.v32.lncRNA_transcripts.fa.gz",
        },
    )

    gencode.data = gencode.data.sample(frac=0.01)
    return gencode


@pytest.fixture
def generate_GENCODE_dask():
    gencode = GENCODE(
        path=
        "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/",
        file_resources={
            "long_noncoding_RNAs.gtf":
            "gencode.v32.long_noncoding_RNAs.gtf.gz",
            "lncRNA_transcripts.fa": "gencode.v32.lncRNA_transcripts.fa.gz",
        },
        npartitions=8,
    )

    gencode.data = gencode.data.sample(frac=0.01)
    return gencode


@pytest.fixture
def generate_MirBase_ftp():
    mirbase = MirBase(path="ftp://mirbase.org/pub/mirbase/CURRENT/")
    mirbase.data = mirbase.data.sample(frac=0.01)
    return mirbase


def test_import_mirbase_db(generate_MirBase_ftp):
    """
    Args:
        generate_MirBase_ftp:
    """
    assert generate_MirBase_ftp.data_path == "ftp://mirbase.org/pub/mirbase/CURRENT/"


def test_import_GENCODE(generate_GENCODE):
    """
    Args:
        generate_GENCODE:
    """
    assert (
        generate_GENCODE.data_path ==
        "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/")


def test_import_dask_GENCODE(generate_GENCODE_dask):
    assert isinstance(generate_GENCODE_dask.data, dd.DataFrame)


def test_annotate_GENCODE(generate_TCGA_LUAD, generate_GENCODE):
    """
    Args:
        generate_TCGA_LUAD:
        generate_GENCODE:
    """
    generate_TCGA_LUAD.LncRNA.annotate_attributes(generate_GENCODE,
                                                  on="gene_id",
                                                  columns=["gene_name"])

    # Test join on off-index
    generate_TCGA_LUAD.LncRNA.annotate_attributes(generate_GENCODE,
                                                  on="gene_name",
                                                  columns=["transcript_id"],
                                                  agg="concat")

    assert {"gene_name", "transcript_id"
            }.issubset(generate_TCGA_LUAD.LncRNA.annotations.columns)


def test_annotate_dask_GENCODE(generate_TCGA_LUAD, generate_GENCODE_dask):
    """
    Args:
        generate_TCGA_LUAD:
        generate_GENCODE_dask:
    """
    generate_GENCODE_dask.data = generate_GENCODE_dask.data[
        generate_GENCODE_dask.data["gene_id"].notnull()]

    # Test join on index column
    generate_TCGA_LUAD.LncRNA.annotate_attributes(generate_GENCODE_dask,
                                                  on="gene_id",
                                                  columns=["gene_name"],
                                                  agg="concat")

    # Test join on off-index
    generate_TCGA_LUAD.LncRNA.annotate_attributes(generate_GENCODE_dask,
                                                  on="gene_name",
                                                  columns=["transcript_id"],
                                                  agg="concat")

    assert {"gene_name", "transcript_id"
            }.issubset(generate_TCGA_LUAD.LncRNA.annotations.columns)


def test_annotate_sequence_GENCODE(generate_TCGA_LUAD, generate_GENCODE):
    generate_TCGA_LUAD.LncRNA.annotate_sequences(generate_GENCODE,
                                                 index="gene_id",
                                                 omic="LncRNA",
                                                 agg="longest")

    assert not generate_TCGA_LUAD.LncRNA.annotations["sequence"].empty
