import io
import os
from glob import glob

import dask.dataframe as dd
import numpy as np
import pandas as pd
# from Bio.UniProt import GOA
from dask import delayed

from .database import Annotatable
from .utils.df import drop_duplicate_columns


class ExpressionData(object):
    def __init__(self,
                 data,
                 transpose,
                 gene_index=None,
                 usecols=None,
                 gene_level=None,
                 sample_level="sample_index",
                 transform_fn=None,
                 dropna=False,
                 npartitions=None,
                 cohort_name=None):
        """
        .. class:: ExpressionData
        This class handles importing of any quantitative omics data that is
        in a table format (e.g. csv, tsv, excel). Pandas will load the DataFrame
        from file with the user-specified columns and genes column name, then
        tranpose it such that the rows are samples and columns are
        gene/transcript/peptides. The user will also specify the index argument,
        which specifies if the genes are ensembl genes ID or gene name, or
        transcripts id/names. The user should be careful about choosing the
        right genes index which makes it easier to annotate functional,
        sequence, and interaction data to it. The dataframe should only contain
        numeric values besides the genes_col_name and the sample barcode id
        indices. :param data: Path or file stream of the table file to import.
        If a pandas DataFrame is passed, then import this dataframe and skip
        preprocessing steps. :type data: str, byte-like, pandas.DataFrame :param
        transpose: True if given data table has samples or columns and variables
        for rows. False if the table has samples for row index, and gene names
        as columns. :type transpose: bool :param gene_index: The column name of
        gene/transcript/protein to index by. :type gene_index: str :param
        usecols: a regex string

            A regex string to import column names from the table. Columns names
            imported are string match, separated by "|".

        Args:
            data:
            transpose:
            gene_index:
            usecols:
            gene_level (str): {"gene", "transcript", "peptide"} Chooses the
                level of the gene/transcript/peptide of the genes list in this
                expression data. The expression DataFrame's index will be
                renamed to this.
            sample_level (str): {"sample_index", "patient_index"} Chooses the
                level of the patient/sample/aliquot indexing.
            transform_fn (bool): default False A callable function to transform
                single values.
            dropna (bool): Whether to drop rows with null values
            npartitions (int): [0-n], default 0 If 0, then uses a Pandas
                DataFrame, if >1, then creates an off-memory Dask DataFrame with
                n partitions
            cohort_name (str): The unique cohort code name string.
        """
        self.cohort_name = cohort_name
        self.gene_level = gene_level
        self.sample_level = sample_level

        df = self.load_dataframe(data,
                                 transpose=transpose,
                                 usecols=usecols,
                                 gene_index=gene_index)
        self.expressions = self.preprocess_table(df,
                                                 usecols=usecols,
                                                 gene_index=gene_index,
                                                 transposed=transpose,
                                                 dropna=dropna)

        # TODO load DD from file directly
        if npartitions and isinstance(self.expressions, pd.DataFrame):
            self.expressions = dd.from_pandas(self.expressions,
                                              npartitions=npartitions)

        if gene_level is not None:
            self.expressions.columns.name = gene_level

        self.expressions.index.name = self.sample_level

        if callable(transform_fn):
            self.expressions = self.expressions.applymap(transform_fn)
        elif transform_fn == "log2":
            self.expressions = self.expressions.applymap(
                lambda x: np.log2(x + 1))

    @property
    def gene_index(self):
        return self.expressions.columns.name

    def load_dataframe(self, data, transpose, usecols, gene_index):
        """Reading table data inputs to create a DataFrame.

        Args:
            data: either a file path, a glob file path (e.g. "table- *.tsv"), a
                pandas.DataFrame, or a dask DataFrame.
            transpose: True if table oriented with samples columns, else False.
            usecols: A regex string to select columns. Default None.
            gene_index:
        """
        if isinstance(data, pd.DataFrame):
            df = data
            df = df.reset_index()
        elif isinstance(data, dd.DataFrame):
            df = data
        elif "*" in data:
            df = self.load_dataframe_glob(data, usecols, gene_index, transpose)
        elif isinstance(data, io.StringIO):
            data.seek(
                0
            )  # Needed since the file was previous read to extract columns information
            df = pd.read_table(data)
        elif type(data) is str and os.path.isfile(data):
            df = pd.read_table(data, sep=None)
        else:
            raise IOError(data)

        # TODO implement handling for multiple file ByteIO
        return df

    def load_dataframe_glob(self, globstring, usecols, genes_index, transpose):
        # type: (str, str, str, bool) -> dd.DataFrame
        """
        Args:
            globstring:
            usecols:
            genes_index:
            transpose:
        """
        lazy_dataframes = []
        for file_path in glob(globstring):
            df = delayed(pd.read_table)(file_path, )
            df = delayed(self.preprocess_table)(df, usecols, genes_index,
                                                transpose, True)
            lazy_dataframes.append(df)

        return dd.from_delayed(lazy_dataframes, divisions=None)

    def preprocess_table(self,
                         df,
                         usecols=None,
                         gene_index=None,
                         transposed=True,
                         sort_index=False,
                         dropna=True):
        # type: (pd.DataFrame, str, str, bool, bool, bool) -> pd.DataFrame
        """This function preprocesses the expression table files where columns
        are samples and rows are gene/transcripts :param df: A Dask or Pandas
        DataFrame :type df: DataFrame :param usecols: A regular expression
        string for the column names to fetch. :type usecols: str :param
        gene_index: The column name containing the gene/transcript names or
        id's. :type gene_index: str :param transposed: Default True. Whether to
        transpose the dataframe so columns are genes (features) and rows are
        samples.

        Args:
            df:
            usecols:
            gene_index:
            transposed:
            sort_index:
            dropna:

        Returns:
            dataframe: a processed Dask DataFrame
        """

        # Filter columns
        if usecols is not None:
            if gene_index not in usecols:
                # include index column in the filter regex query
                usecols = usecols + "|" + gene_index
            df = df.filter(regex=usecols)

        # Drop duplicate sample names
        df = drop_duplicate_columns(df)

        # Drop NA geneID rows
        if dropna:
            df.dropna(axis=0, inplace=True)

        # Remove entries with unknown geneID
        if gene_index is not None:
            df = df[df[gene_index] != '?']
            df.set_index(gene_index, inplace=True)

        # Needed for Dask Delayed
        if sort_index is True:
            df.sort_index(axis=0, ascending=True, inplace=True)

        # Select only numerical columns
        df = df.select_dtypes(include="number")

        # Transpose dataframe to sample rows and gene columns
        if transposed:
            df = df.T

        # Drop duplicate genes
        df = drop_duplicate_columns(df)

        return df

    def set_genes_index(self, index: str, old_index: str):
        """
        Args:
            index (str):
            old_index (str):
        """
        assert isinstance(self, Annotatable) and isinstance(
            self, ExpressionData)
        # Change gene name columns in expressions
        rename_dict = self.get_rename_dict(from_index=old_index,
                                           to_index=index)
        self.expressions.rename(columns=rename_dict, inplace=True)
        self.gene_index = index

        # Change index name in annotation
        self.set_index(index)

    def drop_genes(self, genes_to_drop):
        """
        Args:
            genes_to_drop:
        """
        self.expressions.drop(genes_to_drop, axis=1, inplace=True)
        for gene in genes_to_drop:
            self.features.remove(gene)

    @classmethod
    def name(cls):
        raise NotImplementedError

    def get_genes_list(self, level=None):
        """
        Args:
            level:
        """
        index = self.expressions.columns

        if isinstance(index, pd.MultiIndex):
            return index.get_level_values(
                self.gene_index if level is None else level)
        else:
            return index

    def get_samples_list(self, level=None):
        """
        Args:
            level:
        """
        index = self.expressions.index
        if isinstance(index, pd.MultiIndex):
            return index.get_level_values(
                self.gene_index if level is None else level)
        else:
            return index

    samples = property(get_samples_list)
    features = property(get_genes_list)


class LncRNA(ExpressionData, Annotatable):
    def __init__(self,
                 data,
                 transpose,
                 gene_index=None,
                 usecols=None,
                 gene_level=None,
                 sample_level="sample_index",
                 transform_fn=None,
                 dropna=False,
                 npartitions=None,
                 cohort_name=None):
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
        super(LncRNA, self).__init__(data=data,
                                     transpose=transpose,
                                     gene_index=gene_index,
                                     usecols=usecols,
                                     gene_level=gene_level,
                                     sample_level=sample_level,
                                     transform_fn=transform_fn,
                                     dropna=dropna,
                                     npartitions=npartitions,
                                     cohort_name=cohort_name)

    @classmethod
    def name(cls):
        return cls.__name__


class MessengerRNA(ExpressionData, Annotatable):
    def __init__(self,
                 data,
                 transpose,
                 gene_index=None,
                 usecols=None,
                 gene_level=None,
                 sample_level="sample_index",
                 transform_fn=None,
                 dropna=False,
                 npartitions=None,
                 cohort_name=None):
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
        super(MessengerRNA, self).__init__(data=data,
                                           transpose=transpose,
                                           gene_index=gene_index,
                                           usecols=usecols,
                                           gene_level=gene_level,
                                           sample_level=sample_level,
                                           transform_fn=transform_fn,
                                           dropna=dropna,
                                           npartitions=npartitions,
                                           cohort_name=cohort_name)

    @classmethod
    def name(cls):
        return cls.__name__


class MicroRNA(ExpressionData, Annotatable):
    def __init__(self,
                 data,
                 transpose,
                 gene_index=None,
                 usecols=None,
                 gene_level=None,
                 sample_level="sample_index",
                 transform_fn=None,
                 dropna=False,
                 npartitions=None,
                 cohort_name=None):
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
        super(MicroRNA, self).__init__(data=data,
                                       transpose=transpose,
                                       gene_index=gene_index,
                                       usecols=usecols,
                                       gene_level=gene_level,
                                       sample_level=sample_level,
                                       transform_fn=transform_fn,
                                       dropna=dropna,
                                       npartitions=npartitions,
                                       cohort_name=cohort_name)

    @classmethod
    def name(cls):
        return cls.__name__
