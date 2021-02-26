import copy
import difflib
import gzip
import itertools
import logging
import os
import zipfile
from abc import ABC, abstractmethod

import dask.dataframe as dd
import filetype
import pandas as pd
import rarfile
import validators

from openomics.utils.df import concat_uniques
from openomics.utils.io import get_pkg_data_filename


class Dataset(object):
    COLUMNS_RENAME_DICT = None  # Needs initialization since subclasses may use this field to rename columns in dataframes.

    def __init__(
        self,
        path,
        file_resources=None,
        col_rename=None,
        npartitions=None,
        verbose=False,
    ):
        """This is a class used to instantiate a Dataset given a a set of files from either local files or URLs. When
        creating a Dataset class, the `load_dataframe()` function is called where the file_resources are used to
        load (Pandas or Dask) DataFrames, then performs data wrangling to yield a dataframe at `self.data`. This class
        also provides an interface for -omics tables, e.g. `ExpressionData`, to annotate various annotations, expressions,
        sequences, and disease associations.

        Args:
            path: The folder or url path containing the data file resources. If url path, the files will be downloaded and cached to the user's home folder (at ~/.astropy/).
            file_resources: Used to list required files for preprocessing of the database. A dictionary where keys are required filenames and value are file paths. If None, then the class constructor should automatically build the required file resources dict.
            col_rename (dict): default None, A dictionary to rename columns in the data table. If None, then automatically load defaults.
            npartitions (int): [0-n], default 0 If 0, then uses a Pandas DataFrame, if >1, then creates an off-memory Dask DataFrame with n partitions
            verbose (bool): Default False.
        """
        self.npartitions = npartitions
        self.verbose = verbose

        self.validate_file_resources(path,
                                     file_resources,
                                     npartitions=npartitions,
                                     verbose=verbose)

        self.data = self.load_dataframe(file_resources,
                                        npartitions=npartitions)
        self.data = self.data.reset_index()
        if col_rename is not None:
            self.data = self.data.rename(columns=col_rename)

        self.info() if verbose else None

    def info(self):
        logging.info("{}: {}".format(self.name(), self.data.columns.tolist()))

    def validate_file_resources(self,
                                path,
                                file_resources,
                                npartitions=None,
                                verbose=False) -> None:
        """For each file in file_resources, download the file if path+file is a URL
        or load from disk if a local path. Additionally unzip or unrar if the
        file is compressed.

        Args:
            path (str): The folder or url path containing the data file
                resources. If a url path, the files will be downloaded and cached
                to the user's home folder (at ~/.astropy/).
            file_resources (dict): default None, Used to list required files for
                preprocessing of the database. A dictionary where keys are
                required filenames and value are file paths. If None, then the
                class constructor should automatically build the required file
                resources dict.
            npartitions (int): >0 if the files will be used to create a Dask Dataframe. Default None.
            verbose:
        """
        if validators.url(path):
            for filename, filepath in copy.copy(file_resources).items():
                data_file = get_pkg_data_filename(
                    path, filepath
                )  # Download file and replace the file_resource path
                filepath_ext = filetype.guess(data_file)

                # This null if-clause is needed incase when filetype_ext is None, causing the next clauses to fail
                if filepath_ext is None:
                    file_resources[filename] = data_file

                # Dask will automatically handle uncompression at dd.read_table(compression=filepath_ext)
                elif ".gtf" in filename and npartitions:
                    file_resources[filename] = data_file

                elif filepath_ext.extension == "gz":
                    logging.debug("Uncompressed gzip file at", data_file)
                    file_resources[filename] = gzip.open(data_file, "rt")

                elif filepath_ext.extension == "zip":
                    logging.debug("Uncompressed zip file at", data_file)
                    zf = zipfile.ZipFile(data_file, "r")

                    for subfile in zf.infolist():
                        if (os.path.splitext(subfile.filename)[-1] ==
                                os.path.splitext(filename)[-1]
                            ):  # If the file extension matches
                            file_resources[filename] = zf.open(
                                subfile.filename, mode="r")

                elif filepath_ext.extension == "rar":
                    logging.debug("Uncompressed rar file at", data_file)
                    rf = rarfile.RarFile(data_file, "r")

                    for subfile in rf.infolist():
                        if (os.path.splitext(subfile.filename)[-1] ==
                                os.path.splitext(filename)[-1]
                            ):  # If the file extension matches
                            file_resources[filename] = rf.open(
                                subfile.filename, mode="r")
                else:
                    file_resources[filename] = data_file

        elif os.path.isdir(path) and os.path.exists(path):
            for _, filepath in file_resources.items():
                if not os.path.exists(filepath):
                    raise IOError(filepath)
        else:
            raise IOError(path)

        self.data_path = path
        self.file_resources = file_resources

    def close(self):
        # Close opened file resources
        for filename, filepath in self.file_resources.items():
            if type(self.file_resources[filename]) != str:
                self.file_resources[filename].close()

    @abstractmethod
    def load_dataframe(self, file_resources, npartitions=None):
        # type: (dict, int) -> pd.DataFrame
        """Handles data preprocessing given the file_resources input, and
        returns a DataFrame.

        Args:
            file_resources (dict): A dict with keys as filenames and values as
                full file path.
            npartitions:
        """
        raise NotImplementedError

    @classmethod
    def name(cls):
        return cls.__name__

    @staticmethod
    def list_databases():
        return DEFAULT_LIBRARIES

    def get_annotations(
        self,
        index: str,
        columns: list,
        agg: str = "concat",
        filter_values: pd.Series = None,
    ):
        """Returns the Database's DataFrame such that it's indexed by :param
        index:, which then applies a groupby operation and aggregates all other
        columns by concatenating all unique values.

        Args:
            index (str): The column name of the DataFrame to join by.
            columns ([str]): a list of column names.
            agg (str): Function to aggregate when there is more than one values for each index instance.
                E.g. ['first', 'last', 'sum', 'mean', 'size', 'concat'], default 'concat'.
            filter_values (pd.Series): The values on the `index` column to filter before performing the groupby-agg operations.

        Returns:
            df (DataFrame): A dataframe to be used for annotation
        """
        if not set(columns).issubset(set(self.data.columns)):
            raise Exception(
                "The columns argument must be a list such that it's subset of the following columns in the dataframe",
                "These columns doesn't exist in database:",
                set(columns) - set(self.data.columns.tolist()),
            )

        # Select df columns including df. However the `columns` list shouldn't contain the index column
        if index in columns:
            columns.pop(columns.index(index))

        df = self.data[columns + [index]]

        if filter_values is not None:
            df = df[df[index].isin(list(filter_values))]

        # if index != self.data.index.name and index in self.data.columns:
        #     df = df.set_index(index)

        # Groupby index
        groupby = df.groupby(index)

        #  Aggregate by all columns by concatenating unique values
        if agg == "concat":
            if isinstance(df, pd.DataFrame):
                aggregated = groupby.agg(
                    {col: concat_uniques
                     for col in columns})

            elif isinstance(df, dd.DataFrame):
                collect_concat = dd.Aggregation(
                    name="collect_concat",
                    chunk=lambda s1: s1.apply(list),
                    agg=lambda s2: s2.apply(lambda chunks: filter(
                        lambda x: False if x == "None" or x is None else True,
                        set(itertools.chain.from_iterable(chunks)),
                    )),
                    finalize=lambda s3: s3.apply(lambda xx: "|".join(xx)),
                )
                aggregated = groupby.agg(
                    {col: collect_concat
                     for col in columns})

        # Any other aggregation functions
        else:
            aggregated = groupby.agg({col: agg for col in columns})

        # if aggregated.index.duplicated().sum() > 0:
        #     raise ValueError("DataFrame must not have duplicates in index")
        return aggregated

    def get_expressions(self, index):
        """
        Args:
            index:
        """
        return self.data.groupby(index).median(
        )  # TODO if index by gene, aggregate medians of transcript-level expressions

    @abstractmethod
    def get_rename_dict(self, from_index, to_index):
        """Used to retrieve a lookup dictionary to convert from one index to
        another, e.g., gene_id to gene_name

        Returns
            rename_dict (dict): a rename dict

        Args:
            from_index: an index on the DataFrame for key
            to_index: an index on the DataFrame for value
        """
        raise NotImplementedError


class Annotatable(ABC):
    """This abstract class provides an interface for the omics to annotate external data
    downloaded from various databases. The database will be imported as attributes
    information to the genes's annotations, or interactions between the genes.
    """

    SEQUENCE_COL_NAME = "sequence"

    def get_annotations(self):
        if hasattr(self, "annotations"):
            return self.annotations
        else:
            raise Exception(
                "{} must run initialize_annotations() first.".format(
                    self.name()))

    def get_annotation_expressions(self):
        if hasattr(self, "annotation_expressions"):
            return self.annotation_expressions
        else:
            raise Exception("{} must run annotate_expressions() first.".format(
                self.name()))

    def initialize_annotations(self, index, gene_list=None):
        """
        Args:
            gene_list:
            index:
        """
        if gene_list is None:
            gene_list = self.get_genes_list()

        self.annotations = pd.DataFrame(index=gene_list)
        self.annotations.index.name = index

    def annotate_attributes(self,
                            database: Dataset,
                            on,
                            columns,
                            agg="concat",
                            fuzzy_match=False):
        """Performs a left outer join between the annotation and Database's
        DataFrame, on the index key. The index argument must be column present
        in both DataFrames. If there exists overlapping columns from the join, then
        .fillna() is used to fill NaN values in the old column with non-NaN
        values from the new column.

        Args:
            database (openomics.database.Dataset): Database which contains an dataframe.
            on (str): The column name which exists in both the annotations and Database dataframe to perform the join on.
            columns ([str]): a list of column name to join to the annotation.
            agg (str): Function to aggregate when there is more than one values for each index instance. E.g. ['first', 'last', 'sum', 'mean', 'concat'], default 'concat'.
            fuzzy_match (bool): default False. Whether to join the annotation by applying a fuzzy match on the index with difflib.get_close_matches(). It is very computationally expensive and thus should only be used sparingly.
        """
        if not hasattr(self, "annotations"):
            raise Exception(
                "Must run .initialize_annotations() on, ",
                self.__class__.__name__,
                " first.",
            )

        if on in self.annotations.columns:
            filter_values = self.annotations[on]
        elif on == self.annotations.index.name:
            filter_values = self.annotations.index
        else:
            filter_values = None

        database_df = database.get_annotations(on,
                                               columns=columns,
                                               agg=agg,
                                               filter_values=filter_values)

        if len(database_df.columns) == 0:
            logging.warning(
                "Database annotations is empty and has nothing to annotate.")
            return

        if fuzzy_match:
            database_df.index = database_df.index.map(
                lambda x: difflib.get_close_matches(
                    x, self.annotations.index, n=1)[0])

        # performing join on the index column
        if on == self.annotations.index.name and isinstance(
                database_df, pd.DataFrame):
            new_annotations = self.annotations.join(database_df,
                                                    on=on,
                                                    rsuffix="_")
        # elif on == self.annotations.index.name and isinstance(database_df, dd.DataFrame):
        #     new_annotations = dd.merge(self.annotations, database_df, how="left", on=on, suffixes=("_", ""))

        # performing join on a different column
        else:
            # if isinstance(self.annotations.index, pd.MultiIndex):
            #     old_index = self.annotations.index.names
            # else:
            #     old_index = self.annotations.index.name

            # Save old index, reset the old index, set_index to the join index, perform the join, then change back to the old index
            # This also ensures the index in self.annotations aligns with the gene_index in self.expressions dataframes
            # TODO: Could be very slow on dask dataframes
            # new_annotations = self.annotations.reset_index()
            # new_annotations = new_annotations.set_index(on)
            # new_annotations = new_annotations.join(
            #     database_df, on=on, rsuffix="_").reset_index()
            # new_annotations = new_annotations.set_index(old_index)
            new_annotations = dd.merge(self.annotations,
                                       database_df,
                                       how="left",
                                       on=on,
                                       suffixes=("_", ""))

        # Merge columns if the database DataFrame has overlapping columns with existing column
        duplicate_cols = [
            col for col in new_annotations.columns if col[-1] == "_"
        ]

        # Fill in null values then drop duplicate columns
        for new_col in duplicate_cols:
            old_col = new_col.strip("_")
            new_annotations[old_col].fillna(new_annotations[new_col],
                                            inplace=True,
                                            axis=0)
            new_annotations = new_annotations.drop(columns=new_col)

        # Assign the new results
        self.annotations = new_annotations

    def annotate_sequences(self,
                           database,
                           index,
                           agg="longest",
                           omic=None,
                           **kwargs):
        """Annotate a genes list (based on index) with a dictionary of <gene_name: sequence>. If multiple sequences per
        gene name, then perform some aggregation.

        Args:
            database (Dataset):
            index (str): The gene index column name.
            agg (str): The aggregation method, one of ["longest", "shortest", or "all"]. Default longest.
            omic (str): Default None. Declare the omic type to fetch sequences for.
            **kwargs:
        """
        if omic is None:
            omic = self.name()

        sequences_entries = database.get_sequences(index=index,
                                                   omic=omic,
                                                   agg_sequences=agg,
                                                   **kwargs)

        if type(self.annotations.index) == pd.MultiIndex:
            self.annotations[
                Annotatable.
                SEQUENCE_COL_NAME] = self.annotations.index.get_level_values(
                    index).map(sequences_entries)
        else:
            self.annotations[
                Annotatable.SEQUENCE_COL_NAME] = self.annotations.index.map(
                    sequences_entries)

    def annotate_expressions(self, database, index, fuzzy_match=False):
        """
        Annotate
        Args:
            database:
            index:
            fuzzy_match:
        """
        self.annotation_expressions = pd.DataFrame(
            index=self.annotations.index)

        if self.annotations.index.name == index:
            self.annotation_expressions = self.annotation_expressions.join(
                database.get_expressions(index=index))
        else:
            raise Exception("index argument must be one of",
                            database.data.index)

    def annotate_interactions(self, database, index):
        """
        Args:
            database (Interactions):
            index (str):
        """
        raise NotImplementedError

    def annotate_diseases(self, database, index):
        """
        Args:
            database (DiseaseAssociation):
            index (str):
        """
        self.annotations["disease_associations"] = self.annotations.index.map(
            database.get_disease_assocs(index=index, ))

    def set_index(self, new_index):
        """
        Args:
            new_index:
        """
        self.annotations[new_index].fillna(self.annotations.index.to_series(),
                                           axis=0,
                                           inplace=True)
        self.annotations = self.annotations.reset_index().set_index(new_index)

    def get_rename_dict(self, from_index, to_index):
        """
        Args:
            from_index:
            to_index:
        """
        dataframe = self.annotations.reset_index()
        dataframe = dataframe[dataframe[to_index].notnull()]
        return pd.Series(dataframe[to_index].values,
                         index=dataframe[from_index]).to_dict()


DEFAULT_LIBRARIES = [
    "10KImmunomes"
    "BioGRID"
    "CCLE"
    "DisGeNET"
    "ENSEMBL"
    "GENCODE"
    "GeneMania"
    "GeneOntology"
    "GlobalBiobankEngine"
    "GTEx"
    "HMDD_miRNAdisease"
    "HPRD_PPI"
    "HUGO_Gene_names"
    "HumanBodyMapLincRNAs"
    "IntAct"
    "lncBase"
    "LNCipedia"
    "LncReg"
    "lncRInter"
    "lncrna2target"
    "lncRNA_data_repository"
    "lncrnadisease"
    "lncRNome"
    "mirbase"
    "miRTarBase"
    "NHLBI_Exome_Sequencing_Project"
    "NONCODE"
    "NPInter"
    "PIRD"
    "RegNetwork"
    "RISE_RNA_Interactions"
    "RNAcentral"
    "StarBase_v2.0"
    "STRING_PPI"
    "TargetScan"
]
