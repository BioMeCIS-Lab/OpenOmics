import networkx as nx
import pandas as pd

from .database import Annotatable
from .transcriptomics import Expression


class Protein(Expression, Annotatable):
    def __init__(
        self,
        data,
        transpose,
        gene_index=None,
        usecols=None,
        gene_level=None,
        sample_level="sample_index",
        transform_fn=None,
        dropna=False,
        npartitions=None,
        cohort_name=None,
    ):
        """
        Args:
            data:
            transpose:
            gene_index:
            usecols:
            gene_level:
            sample_level:
            transform_fn:
            dropna:
            npartitions:
            cohort_name:
        """
        super(Protein, self).__init__(
            data=data,
            transpose=transpose,
            gene_index=gene_index,
            usecols=usecols,
            gene_level=gene_level,
            sample_level=sample_level,
            transform_fn=transform_fn,
            dropna=dropna,
            npartitions=npartitions,
            cohort_name=cohort_name,
        )

    @classmethod
    def name(cls):
        return cls.__name__

    def process_HPRD_PPI_network(self, ppi_data_file_path):
        """
        Args:
            ppi_data_file_path:
        """
        HPRD_PPI = pd.read_table(ppi_data_file_path, header=None)
        self.HPRD_PPI_network = nx.from_pandas_edgelist(
            HPRD_PPI, source=0, target=3, create_using=nx.DiGraph())
