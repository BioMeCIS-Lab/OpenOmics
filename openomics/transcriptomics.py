from typing import Union, Dict
import logging
import io
import re
import os
import validators
from glob import glob
from typing import Union

import dask.dataframe as dd
import numpy as np
import pandas as pd

# from Bio.UniProt import GOA
from dask import delayed

from .database import Annotatable
from .utils.df import drop_duplicate_columns
from .utils import get_pkg_data_filename


class Expression(object):
    """This class handles importing of any quantitative omics data that is
    in a table format (e.g. csv, tsv, excel). Pandas will load the DataFrame
    from file with the user-specified columns and genes column name, then
    tranpose it such that the rows are samples and columns are
    gene/transcript/peptides. The user will also specify the index argument,
    which specifies if the genes are ensembl genes ID or gene name, or
    transcripts id/names. The user should be careful about choosing the
    right genes index which makes it easier to annotate functional,
    sequence, and interaction data to it. The dataframe should only contain
    numeric values besides the genes_col_name and the sample barcode id
    indices.

    """
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
        prefix=None,
        **kwargs,
    ):
        """This class handles importing of any quantitative omics data that is
        in a table format (e.g. csv, tsv, excel). Pandas will load the DataFrame
        from file with the user-specified columns and genes column name, then
        tranpose it such that the rows are samples and columns are
        gene/transcript/peptides. The user will also specify the index argument,
        which specifies if the genes are ensembl genes ID or gene name, or
        transcripts id/names. The user should be careful about choosing the
        right genes index which makes it easier to annotate functional,
        sequence, and interaction data to it. The dataframe should only contain
        numeric values besides the genes_col_name and the sample barcode id
        indices.

        Args:
            data (str, byte-like, pandas.DataFrame): Path or file stream of the
                table file to import. If a pandas DataFrame is passed, then
                import this dataframe and skip preprocessing steps.
            transpose (bool): True if given data table has samples or columns
                and variables for rows. False if the table has samples for row
                index, and gene names as columns.
            gene_index (str): The column name of gene/transcript/protein to
                index by.
            usecols: A regex string to import column names from the table.
                Columns names imported are string match, separated by "|".
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
        """
        self.gene_level = gene_level
        self.sample_level = sample_level

        df = self.load_dataframe(
            data,
            transpose=transpose,
            usecols=usecols,
            gene_index=gene_index,
            dropna=dropna,
        )
        self.expressions = self.preprocess_table(
            df,
            usecols=usecols,
            gene_index=gene_index,
            transposed=transpose,
            dropna=dropna,
        )

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

    def load_dataframe(
        self,
        data: Union[str, pd.DataFrame, dd.DataFrame, io.StringIO],
        transpose: bool,
        usecols: str,
        gene_index: str,
        dropna: bool,
    ):
        """Reading table data inputs to create a DataFrame.

        Args:
            data: either a file path, a glob file path (e.g. "table-*.tsv"), a
                pandas.DataFrame, or a dask DataFrame.
            transpose (bool): True if table oriented with samples columns, else
                False.
            usecols (str): A regex string to select columns. Default None.
            gene_index (str):
            dropna (bool): Whether to drop rows with null values

        Returns:
            Union[pd.DataFrame, dd.DataFrame]: The preprocessed dataframe.
        """
        if isinstance(data, pd.DataFrame):
            df = data

        elif isinstance(data, dd.DataFrame):
            df = data

        elif isinstance(data, str) and "*" in data:
            # TODO implement handling for multiple file ByteIO
            df = self.load_dataframe_glob(
                globstring=data,
                usecols=usecols,
                gene_index=gene_index,
                transpose=transpose,
                dropna=dropna,
            )

        elif isinstance(data, io.StringIO):
            # Needed since the file was previous read to extract columns information
            data.seek(0)
            df = pd.read_table(data)

        elif isinstance(data, str) and os.path.isfile(data):
            df = pd.read_table(data, sep=None, engine="python")

        elif isinstance(data, str) and validators.url(data):
            dataurl, filename = os.path.split(data)
            file = get_pkg_data_filename(dataurl + "/", filename)
            df = pd.read_table(file)

        else:
            raise IOError(data)

        return df

    def preprocess_table(
        self,
        df: Union[pd.DataFrame, dd.DataFrame],
        usecols: str = None,
        gene_index: str = None,
        transposed: bool = True,
        sort_index: bool = False,
        dropna: bool = True,
    ):
        """This function preprocesses the expression table files where columns
        are samples and rows are gene/transcripts :param df: A Dask or Pandas
        DataFrame :type df: DataFrame :param usecols: A regular expression
        string for the column names to fetch. :type usecols: str :param
        gene_index: The column name containing the gene/transcript names or
        id's. :type gene_index: str :param transposed: Default True. Whether to
        transpose the dataframe so columns are genes (features) and rows are
        samples.

        Args:
            df (pd.DataFrame):
            usecols (str):
            gene_index (str):
            transposed (bool):
            sort_index (bool):
            dropna (bool):
        Returns:
            dataframe: a processed Dask DataFrame
        """
        # Filter columns
        if usecols is not None and isinstance(usecols, str):
            if gene_index not in usecols:
                # include index column in the filter regex query
                usecols = usecols + "|" + gene_index

            if isinstance(df, pd.DataFrame):
                df = df.filter(regex=usecols)
            elif isinstance(df, dd.DataFrame):
                r = re.compile(usecols)
                columns = list(filter(r.match, df.columns))
                df = df[columns]

        elif usecols is not None and isinstance(usecols, list):
            if gene_index not in usecols:
                usecols.append(gene_index)
            df = df[usecols]

        # Drop duplicate sample names
        df = drop_duplicate_columns(df)

        # Drop NA geneID rows
        if dropna:
            df.dropna(axis=0, inplace=True)

        if gene_index is not None and df.index.name != gene_index:
            print(gene_index, "index", df.index.name, "cols", df.columns)
            df = df.set_index(gene_index)

        # Needed for Dask Delayed
        if sort_index is True:
            df = df.sort_index(axis=0, ascending=True)

        # Select only numerical columns
        df = df.select_dtypes(include="number")

        # Transpose dataframe to sample rows and gene columns
        if transposed:
            df = df.T

        # Drop duplicate genes
        df = drop_duplicate_columns(df)

        return df

    def load_dataframe_glob(
        self,
        globstring: str,
        usecols: str,
        gene_index: str,
        transpose: bool,
        dropna: bool,
    ):
        """
        Args:
            globstring (str):
            usecols (str):
            gene_index (str):
            transpose (bool):
        Returns:
            dd.DataFrame
        """
        def convert_numerical_to_float(df: pd.DataFrame):
            cols = df.columns[~df.dtypes.eq("object")]
            df[cols] = df[cols].astype(float)
            return df

        filenames = []

        lazy_dataframes = []
        for file_path in glob(globstring):
            filenames.append(os.path.split(file_path)[1])

            df = delayed(pd.read_table)(file_path)
            # df = delayed(convert_numerical_to_float)(df)
            df = delayed(self.preprocess_table)(
                df,
                usecols,
                gene_index,
                transpose,
                True,
                dropna  # sort_index
            )
            lazy_dataframes.append(df)

        logging.info("Files matched: {}".format(filenames))

        return dd.from_delayed(lazy_dataframes,
                               divisions=None,
                               verify_meta=True)

    def set_genes_index(self, index: str, old_index: str):
        """
        Args:
            index (str):
            old_index (str):
        """
        assert isinstance(self, Annotatable) and isinstance(self, Expression)
        # Change gene name columns in expressions
        rename_dict = self.get_rename_dict(from_index=old_index,
                                           to_index=index)
        self.expressions.rename(columns=rename_dict, inplace=True)
        self.gene_index = index

        # Change index name in annotation
        self.set_index(index)

    def drop_genes(self, gene_ids: str):
        """Drop columns representing genes/rna/proteins in self.expressions
        dataframe.

        Args:
            gene_ids (str): list of strings that are a subset of the columns
                list
        """
        self.expressions = self.expressions.drop(gene_ids, axis=1)
        if hasattr(self, "annotations") and not self.annotations.empty:
            self.annotations = self.annotations.drop(gene_ids, axis=0)

    def drop_samples(self, sample_ids):
        """
        Args:
            sample_ids:
        """
        self.expressions = self.expressions.drop(sample_ids, axis=0)

    @classmethod
    def name(cls):
        raise NotImplementedError

    def get_genes_list(self, level: int = None):
        """
        Args:
            level (int): Default None. Only needed if gene index is a :class:`pd.MultiIndex`
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


class LncRNA(Expression, Annotatable):
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
        super(LncRNA, self).__init__(
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


class MessengerRNA(Expression, Annotatable):
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
        super(MessengerRNA, self).__init__(
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


class MicroRNA(Expression, Annotatable):
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
        super(MicroRNA, self).__init__(
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
