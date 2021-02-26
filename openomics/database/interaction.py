from abc import abstractmethod
import logging

import networkx as nx
from Bio import SeqIO

from openomics.database.annotation import *
from openomics.database.base import Dataset
from openomics.database.sequence import SequenceDataset


class Interactions(Dataset):
    def __init__(
        self,
        path,
        file_resources,
        source_col_name=None,
        target_col_name=None,
        source_index=None,
        target_index=None,
        edge_attr=None,
        filters=None,
        directed=True,
        relabel_nodes=None,
        verbose=False,
    ):
        """
        This is an abstract class used to instantiate a database given a folder containing various file resources. When creating a Database class, the load_data function is called where the file resources are load as a DataFrame and performs necessary processings. This class provides an interface for RNA classes to annotate various genomic annotation, functional annotation, sequences, and disease associations.
        Args:
            path (str):
                The folder path containing the data files.
            file_resources (dict): default None,
                Used to list required files for load_network of the dataset. A dictionary where keys are required filenames and value are file paths. If None, then the class constructor should automatically build the required file resources dict.
            source_col_name (str):
                Column name of DataFrame to be used as the source node names.
            target_col_name (str):
                Column name of DataFrame to be used as the target node names.
            source_index (str):
                One of {"gene_name", "gene_id", "transcript_name", "transcript_id", "protein_name", "protein_id"}
            target_index (str):
                One of {"gene_name", "gene_id", "transcript_name", "transcript_id", "protein_name", "protein_id"}
            edge_attr (list):
                A list of column names to be included as attributes for each edge (source-target pairs).
            filters (dict):
                Optional. A dict with key matching the data table (from load_network()) columns and values for the filtering on that column.
            directed (bool): default True,
                Whether to create a directed or an undirected network.
            relabel_nodes (dict): default None,
                A dictionary to rename nodes in the network, where the nodes with name <dict[key]> will be renamed to <dict[value]>
        """
        # This class should NOT call super's __init__()
        self.validate_file_resources(path, file_resources, verbose=verbose)

        self.data_path = path
        self.file_resources = file_resources
        self.source_index = source_index
        self.target_index = target_index
        self.network = self.load_network(
            file_resources=file_resources,
            source_col_name=source_col_name,
            target_col_name=target_col_name,
            edge_attr=edge_attr,
            directed=directed,
            filters=filters,
        )
        assert isinstance(self.network, nx.Graph), (
            f"load_network must return a NetworkX Graph or DiGraph but instead "
            f"returns {self.network}")
        self.network.name = self.name()

        if self.network is None:
            raise Exception(
                "Make sure load_network() returns a Networkx Graph and is called with super().__init__() in the constructor."
            )

        if relabel_nodes is not None:
            self.network = nx.relabel_nodes(self.network,
                                            mapping=relabel_nodes)

        self.verbose = verbose
        self.info() if verbose else None

    def info(self):
        print("{}".format(nx.info(self.network)))

    @classmethod
    def name(cls):
        return cls.__name__

    @abstractmethod
    def load_network(
        self,
        file_resources,
        source_col_name,
        target_col_name,
        edge_attr,
        directed,
        filters,
    ):
        """
        Handles data processing from `file_resources` to a Pandas DataFrame which contain edgelist data, then constructs
        and return a NetworkX Graph.
        Args:
            file_resources: a dict of file name and file path/object
            source_col_name (str): column name of the dataframe for source in the edge
            target_col_name (str): column name of the dataframe for target in the edge
            edge_attr (list): list of str for column data to include in each edge
            directed (bool): True to return a DiGraph(), else Graph()
            filters: A dict of {column name: column values} to filter the dataframe
        Returns:
            network: a NetworkX Graph or DiGraph
        """
        raise NotImplementedError

    def get_interactions(self, nodelist=None, data=False, inclusive=True):
        """

        Args:
            nodelist (list):
                A list of nodes to fetch edges from
            data (bool): default False
                Whether to include edge attributes
            inclusive (bool): default False
                Whether to only retrieve edges from nodes inclusive in nodelist.

        Returns:
            edges (OutEdgeView): a NetworkX edgelist
        """
        if not hasattr(self, "network"):
            raise Exception(
                "{} does not have network interaction data yet. Must run load_network() and assign self.network field first."
                .format(self.name()))

        if nodelist is None:
            return self.network.edges(data=data)

        if inclusive:
            return self.network.subgraph(nodelist).edges(data=data)
        else:
            return self.network.edges(nbunch=nodelist, data=data)

    def filter_values(self, df: pd.DataFrame, filters: dict, case=False):
        if filters is None:
            return df

        for key, values in filters.items():
            if key not in df.columns:
                logging.info(
                    f"Filter key `{key}` must be in one of {df.columns}")
                continue
            n_rows = df.shape[0]

            if isinstance(values, list):
                if case is False:
                    df = df[df[key].str.upper().isin(
                        [val.upper() for val in values])]
                else:
                    df = df[df[key].isin(values)]
            elif isinstance(values, str):
                df = df[df[key].str.contains(values, case=case)]
            else:
                df = df[df[key] == values]

            logging.info(
                f"INFO: Removed {n_rows - df.shape[0]} rows with `{key}` != {values}"
            )

        assert (df.shape[0] >
                0), f"ERROR: Dataframe is empty because of filter: {filters}"
        return df


class GeneMania(Interactions):
    def __init__(
        self,
        path,
        file_resources=None,
        source_col_name="Gene_A",
        target_col_name="Gene_B",
        source_index="gene_name",
        target_index="gene_name",
        edge_attr=None,
        filters=None,
        directed=True,
        relabel_nodes=None,
    ):
        if edge_attr is None:
            edge_attr = ["Weight"]
        if file_resources is None:
            file_resources = {}
            file_resources[
                "COMBINED.DEFAULT_NETWORKS.BP_COMBINING.txt"] = os.path.join(
                    path, "COMBINED.DEFAULT_NETWORKS.BP_COMBINING.txt")
            file_resources["identifier_mappings.txt"] = os.path.join(
                path, "identifier_mappings.txt")

        super(GeneMania, self).__init__(
            path=path,
            file_resources=file_resources,
            source_col_name=source_col_name,
            target_col_name=target_col_name,
            source_index=source_index,
            target_index=target_index,
            edge_attr=edge_attr,
            filters=filters,
            directed=directed,
            relabel_nodes=relabel_nodes,
        )

    def load_network(
        self,
        file_resources,
        source_col_name,
        target_col_name,
        edge_attr,
        directed,
        filters,
    ):
        interactions = pd.read_table(
            file_resources["COMBINED.DEFAULT_NETWORKS.BP_COMBINING.txt"],
            low_memory=True,
        )
        identifier = pd.read_table(file_resources["identifier_mappings.txt"])

        # Rename ENSG ID's to gene names
        identifier = identifier[identifier["Source"] == "Gene Name"]
        id_mapping = pd.Series(identifier["Name"].values,
                               index=identifier["Preferred_Name"]).to_dict()
        interactions.replace(id_mapping, inplace=True)

        genemania_RNA_RNA_network = nx.from_pandas_edgelist(
            interactions,
            source=source_col_name,
            target=target_col_name,
            edge_attr=edge_attr,
            create_using=nx.DiGraph(),
        )
        return genemania_RNA_RNA_network


class BioGRID(Interactions):
    def __init__(
        self,
        path="https://downloads.thebiogrid.org/Download/BioGRID/Latest-Release/",
        file_resources=None,
        source_col_name="Official Symbol Interactor A",
        target_col_name="Official Symbol Interactor B",
        source_index="gene_name",
        target_index="gene_name",
        edge_attr=[
            "Score",
            "Throughput",
            "Experimental System",
            "Experimental System Type",
        ],
        filters={"Organism Interactor A": 9606},
        directed=False,
        relabel_nodes=None,
    ):
        if file_resources is None:
            file_resources = {}
            file_resources["BIOGRID-ALL-LATEST.tab2.zip"] = os.path.join(
                path, "BIOGRID-ALL-LATEST.tab2.zip")

        super(BioGRID, self).__init__(
            path=path,
            file_resources=file_resources,
            source_col_name=source_col_name,
            target_col_name=target_col_name,
            source_index=source_index,
            target_index=target_index,
            edge_attr=edge_attr,
            directed=directed,
            relabel_nodes=relabel_nodes,
            filters=filters,
        )

    def load_network(
        self,
        file_resources,
        source_col_name,
        target_col_name,
        edge_attr,
        directed,
        filters,
    ):
        df = pd.read_table(
            file_resources["BIOGRID-ALL-LATEST.tab2.zip"],
            na_values=["-"],
            # usecols=['Official Symbol Interactor A', 'Official Symbol Interactor B',
            #          'Organism Interactor A', 'Score', 'Throughput', 'Qualifications',
            #          'Modification', 'Phenotypes'],
            low_memory=True,
        )

        logging.info("{}: {}".format(self.name(), df.columns.tolist()))

        df = self.filter_values(df, filters)
        network = nx.from_pandas_edgelist(
            df,
            source=source_col_name,
            target=target_col_name,
            edge_attr=edge_attr,
            create_using=nx.DiGraph() if directed else nx.Graph(),
        )
        return network


class STRING(Interactions, SequenceDataset):
    COLUMNS_RENAME_DICT = {
        "protein_external_id": "protein_id",
        "preferred_name": "protein_name",
    }

    def __init__(
        self,
        path="https://stringdb-static.org/download/",
        file_resources=None,
        species_id="9606",
        source_col_name="item_id_a",
        target_col_name="item_id_b",
        source_index="protein_name",
        target_index="protein_name",
        edge_attr=["score"],
        directed=False,
        relabel_nodes=None,
        verbose=False,
    ):
        """

        Args:
            species_id (str): Required. Must provide species id number to download the correct STRING dataset.
        """
        if file_resources is None:
            file_resources = {}
            file_resources["protein.actions.txt"] = os.path.join(
                path,
                "protein.actions.v11.0/{}.protein.actions.v11.0.txt.gz".format(
                    species_id),
            )
            file_resources["protein.links.txt"] = os.path.join(
                path,
                "protein.links.v11.0/{}.protein.links.v11.0.txt.gz".format(
                    species_id),
            )
            file_resources["protein.info.txt"] = os.path.join(
                path,
                "protein.info.v11.0/{}.protein.info.v11.0.txt.gz".format(
                    species_id),
            )
            file_resources["protein.sequences.fa"] = os.path.join(
                path,
                "protein.sequences.v11.0/{}.protein.sequences.v11.0.fa.gz".
                format(species_id),
            )

        super(STRING, self).__init__(
            path=path,
            file_resources=file_resources,
            source_col_name=source_col_name,
            target_col_name=target_col_name,
            source_index=source_index,
            target_index=target_index,
            edge_attr=edge_attr,
            directed=directed,
            relabel_nodes=relabel_nodes,
            verbose=verbose,
        )

        self.file_resources["protein.info.txt"].seek(0)
        self.df = pd.read_table(file_resources["protein.info.txt"])
        self.df = self.df.reset_index()
        self.df = self.df.rename(columns=self.COLUMNS_RENAME_DICT)

    def load_network(
        self,
        file_resources,
        source_col_name,
        target_col_name,
        edge_attr,
        directed,
        filters,
    ):
        # protein_interactions = pd.read_table(file_resources["protein.links.txt"], sep=" ", low_memory=True)
        protein_interactions = pd.read_table(
            file_resources["protein.actions.txt"], sep="\t", low_memory=True)
        logging.info("{}: {}".format(self.name(),
                                     protein_interactions.columns.tolist()))
        protein_info = pd.read_table(file_resources["protein.info.txt"])

        self.protein_id2name = protein_info.set_index(
            "protein_external_id")["preferred_name"].to_dict()
        network = nx.from_pandas_edgelist(
            protein_interactions,
            source=source_col_name,
            target=target_col_name,
            edge_attr=edge_attr,
            create_using=nx.DiGraph() if directed else nx.Graph(),
        )
        network = nx.relabel_nodes(network, self.protein_id2name)
        return network

    def get_sequences(self,
                      index="protein_name",
                      omic=None,
                      agg_sequences=None):
        if hasattr(self, "seq_dict"):
            return self.seq_dict

        self.seq_dict = {}
        collisions = 0
        for record in SeqIO.parse(self.file_resources["protein.sequences.fa"],
                                  "fasta"):
            gene_id = str(record.name)
            gene_name = self.protein_id2name[gene_id]
            sequence_str = str(record.seq)
            if index == "protein_name":
                key = gene_name
            elif index == "protein_id":
                key = gene_id

            if key in self.seq_dict:
                collisions += 1

            self.seq_dict[key] = sequence_str

        logging.info("Seq {} collisions: {}".format(index, collisions))
        return self.seq_dict


class LncBase(Interactions, Dataset):
    def __init__(
        self,
        path,
        file_resources=None,
        strip_mirna_name=False,
        source_col_name="mirna",
        target_col_name="geneId",
        source_index="transcript_name",
        target_index="gene_id",
        edge_attr=None,
        filters={"species": "Homo sapiens"},
        directed=True,
        relabel_nodes=None,
    ):
        self.strip_mirna_name = strip_mirna_name

        if edge_attr is None:
            edge_attr = ["tissue", "positive_negative"]
        if file_resources is None:
            file_resources = {}
            file_resources["LncBasev2_download.csv"] = os.path.join(
                path, "LncBasev2_download.csv")

        super(LncBase, self).__init__(
            path=path,
            file_resources=file_resources,
            source_col_name=source_col_name,
            target_col_name=target_col_name,
            source_index=source_index,
            target_index=target_index,
            edge_attr=edge_attr,
            directed=directed,
            relabel_nodes=relabel_nodes,
            filters=filters,
        )

    def get_rename_dict(self, from_index="geneId", to_index="geneName"):
        lncbase_df = pd.read_table(
            self.file_resources["LncBasev2_download.csv"], low_memory=True)
        gene_id_to_gene_name_dict = pd.Series(
            lncbase_df["geneName"].values,
            index=lncbase_df["geneId"]).to_dict()
        return gene_id_to_gene_name_dict

    def load_network(
        self,
        file_resources,
        source_col_name="mirna",
        target_col_name="gene_id",
        edge_attr=None,
        directed=True,
        filters=None,
    ):
        if edge_attr is None:
            edge_attr = ["tissue", "positive_negative"]
        df = pd.read_table(file_resources["LncBasev2_download.csv"],
                           low_memory=True)
        logging.info(f"{self.name()}, {df.columns.tolist()}")
        df.replace(
            {
                "species": {
                    "Homo Sapiens": "Homo sapiens",
                    "Mus Musculus": "Mus musculus",
                }
            },
            inplace=True,
        )

        df = self.filter_values(df, filters)

        if self.strip_mirna_name:
            df["mirna"] = df["mirna"].str.lower()
            df["mirna"] = df["mirna"].str.replace("-3p.*|-5p.*", "")

        lncBase_lncRNA_miRNA_network = nx.from_pandas_edgelist(
            df,
            source=source_col_name,
            target=target_col_name,
            edge_attr=edge_attr,
            create_using=nx.DiGraph() if directed else nx.Graph(),
        )
        return lncBase_lncRNA_miRNA_network


class LncReg(Interactions):
    def __init__(
        self,
        path,
        file_resources,
        source_col_name="A_name_in_paper",
        target_col_name="B_name_in_paper",
        source_index="transcript_name",
        target_index="gene_name",
        edge_attr=["relationship", "mechanism", "pmid"],
        filters=None,
        directed=True,
        relabel_nodes=None,
        verbose=False,
    ):
        if file_resources is None:
            file_resources = {}
            file_resources["data.xlsx"] = os.path.join(path, "data.xlsx")

        super(LncReg, self).__init__(
            path,
            file_resources=file_resources,
            source_col_name=source_col_name,
            target_col_name=target_col_name,
            source_index=source_index,
            target_index=target_index,
            edge_attr=edge_attr,
            filters=filters,
            directed=directed,
            relabel_nodes=relabel_nodes,
            verbose=verbose,
        )

    def load_network(
        self,
        file_resources,
        source_col_name,
        target_col_name,
        edge_attr,
        directed,
        filters,
    ):
        df = pd.read_excel(self.file_resources["data.xlsx"])
        logging.info(f"{self.name()}, {df.columns.tolist()}")

        df = df[df["species"] == "Homo sapiens"]
        df.loc[df["B_category"] == "miRNA", "B_name_in_paper"] = df[
            df["B_category"] == "miRNA"]["B_name_in_paper"].str.replace(
                "-3p.*|-5p.*", "")
        df.loc[df["B_category"] == "miRNA", "B_name_in_paper"] = df[
            df["B_category"] == "miRNA"]["B_name_in_paper"].str.replace(
                "MIR", "hsa-mir-")
        df.loc[df["B_category"] == "miRNA", "B_name_in_paper"] = df[
            df["B_category"] == "miRNA"]["B_name_in_paper"].str.replace(
                "let-", "hsa-let-")

        LncReg_lncRNA_RNA_network = nx.from_pandas_edgelist(
            df,
            source=source_col_name,
            target=target_col_name,
            edge_attr=edge_attr,
            create_using=nx.DiGraph(),
        )
        return LncReg_lncRNA_RNA_network


class lncRInter(Interactions):
    def __init__(
        self,
        path,
        file_resources=None,
        source_col_name="lncrna",
        target_col_name="Interacting partner",
        source_index="gene_name",
        target_index="gene_name",
        edge_attr=None,
        filters=None,
        directed=True,
        relabel_nodes=None,
    ):
        if edge_attr is None:
            edge_attr = [
                "Interaction Class", "Interaction Mode", "Tissue", "Phenotype"
            ]
        if file_resources is None:
            file_resources = {}
            file_resources["human_interactions.txt"] = os.path.join(
                path, "human_interactions.txt")

        super(lncRInter, self).__init__(
            path,
            file_resources,
            source_col_name=source_col_name,
            target_col_name=target_col_name,
            source_index=source_index,
            target_index=target_index,
            edge_attr=edge_attr,
            directed=directed,
            relabel_nodes=relabel_nodes,
            filters=filters,
        )

    def load_network(
        self,
        file_resources,
        source_col_name,
        target_col_name,
        edge_attr,
        directed,
        filters,
    ):
        lncRInter_df = pd.read_table(file_resources["human_interactions.txt"])
        logging.info(f"{self.name()}, {lncRInter_df.columns.tolist()}")

        lncRInter_df = self.filter_values(lncRInter_df, filters)
        # Data cleaning
        lncRInter_df.loc[
            lncRInter_df["Interacting partner"].str.contains("MIR"),
            "Interacting partner", ] = lncRInter_df.loc[
                lncRInter_df["Interacting partner"].str.contains("MIR"),
                "Interacting partner", ].str.lower()
        lncRInter_df["Interacting partner"] = lncRInter_df[
            "Interacting partner"].str.replace("mirlet", "hsa-let-")
        lncRInter_df["Interacting partner"] = lncRInter_df[
            "Interacting partner"].str.replace("mir", "hsa-mir-")
        lncRInter_df["Interacting partner"][
            lncRInter_df["Interacting partner"].str.contains(
                r"[mir|let]\-[\d]+[a-z]+[\d]+"
            )] = lncRInter_df["Interacting partner"][
                lncRInter_df["Interacting partner"].str.contains(
                    r"[mir|let]\-[\d]+[a-z]+[\d]+")].apply(
                        lambda x: x[:-1] + "-" + x[-1])

        lncRInter_network = nx.from_pandas_edgelist(
            lncRInter_df,
            source=source_col_name,
            target=target_col_name,
            edge_attr=edge_attr,
            create_using=nx.DiGraph() if directed else nx.Graph(),
        )
        return lncRInter_network


class LncRNA2Target(Interactions):
    def __init__(
        self,
        path="http://123.59.132.21/lncrna2target/data/",
        file_resources=None,
        source_index="gene_name",
        target_index="gene_name",
        edge_attr=None,
        filters={
            "species_id": 9606,
            "Species": "Homo sapiens"
        },
        directed=True,
        relabel_nodes=None,
        version="high_throughput",
    ):
        """

        Args:
            version (str): one of ["high_throughput", "low_throughput"].
                The high_throughput version of lncRNA2Target database is v2.0 and low_throughput is v1.0, according to the database's website.
            species_id (str, int): one of [9606, "Homo sapiens"].
                The species column in high_throughput is formatted in int (e.g. 9606) and in low_throughput is in str (e.g. "Homo sapiens")
        """
        self.version = version
        if file_resources is None:
            file_resources = {}
            file_resources[
                "lncRNA_target_from_high_throughput_experiments.txt"] = os.path.join(
                    path, "lncrna_target.rar")
            file_resources[
                "lncRNA_target_from_low_throughput_experiments.xlsx"] = os.path.join(
                    path, "lncRNA_target_from_low_throughput_experiments.xlsx")

        if self.version == "high_throughput":
            super(LncRNA2Target, self).__init__(
                path,
                file_resources,
                source_col_name="lncrna_symbol",
                target_col_name="gene_symbol",
                source_index=source_index,
                target_index=target_index,
                edge_attr=edge_attr,
                filters=filters,
                directed=directed,
                relabel_nodes=relabel_nodes,
            )
        if self.version == "low_throughput":
            super(LncRNA2Target, self).__init__(
                path,
                file_resources,
                source_col_name="GENCODE_gene_name",
                target_col_name="Target_official_symbol",
                source_index=source_index,
                target_index=target_index,
                edge_attr=edge_attr,
                filters=filters,
                directed=directed,
                relabel_nodes=relabel_nodes,
            )

    def load_network(
        self,
        file_resources,
        source_col_name,
        target_col_name,
        edge_attr,
        directed,
        filters,
    ):
        if self.version == "high_throughput":
            return self.load_network_high_throughput(file_resources,
                                                     source_col_name,
                                                     target_col_name,
                                                     edge_attr, directed)
        elif self.version == "low_throughput":
            return self.load_network_low_throughput(file_resources,
                                                    source_col_name,
                                                    target_col_name, edge_attr,
                                                    directed)
        else:
            raise Exception(
                "LncRNA2Target version argument must be one of 'high_throughput' or 'low_throughput'"
            )

    def load_network_high_throughput(
        self,
        file_resources,
        source_col_name="lncrna_symbol",
        target_col_name="gene_symbol",
        edge_attr=None,
        directed=True,
        filters=None,
    ):
        table = pd.read_table(
            file_resources[
                "lncRNA_target_from_high_throughput_experiments.txt"],
            sep="\t",
        )
        table = self.filter_values(table, filters)
        logging.info(f"{self.name()}, {table.columns.tolist()}")

        table["lncrna_symbol"] = table["lncrna_symbol"].str.upper()
        table["lncrna_symbol"] = table["lncrna_symbol"].str.replace("LINC", "")
        table["gene_symbol"] = table["gene_symbol"].str.upper()
        lncrna2target_high_throughput_network = nx.from_pandas_edgelist(
            table,
            source=source_col_name,
            target=target_col_name,
            edge_attr=edge_attr,
            create_using=nx.DiGraph() if directed else nx.Graph(),
        )
        return lncrna2target_high_throughput_network

    def load_network_low_throughput(
        self,
        file_resources,
        source_col_name="GENCODE_gene_name",
        target_col_name="Target_official_symbol",
        edge_attr=None,
        directed=True,
        filters=None,
    ):
        table = pd.read_excel(file_resources[
            "lncRNA_target_from_low_throughput_experiments.xlsx"])
        table = self.filter_values(table, filters)
        logging.info(f"{self.name()}, {table.columns.tolist()}")

        table["Target_official_symbol"] = table[
            "Target_official_symbol"].str.replace("(?i)(mir)",
                                                  "hsa-mir-",
                                                  regex=True)
        table["Target_official_symbol"] = table[
            "Target_official_symbol"].str.replace("--", "-")
        table["Target_official_symbol"].apply(
            lambda x: x.lower() if "mir" in x.lower() else x.upper())
        table["GENCODE_gene_name"] = table["GENCODE_gene_name"].str.upper()
        lncrna2target_low_throughput_network = nx.from_pandas_edgelist(
            table,
            source=source_col_name,
            target=target_col_name,
            edge_attr=edge_attr,
            create_using=nx.DiGraph() if directed else nx.Graph(),
        )
        return lncrna2target_low_throughput_network


class lncRNome(Interactions, Dataset):
    def __init__(
        self,
        path,
        file_resources,
        source_col_name="Gene Name",
        target_col_name="Binding miRNAs",
        source_index="gene_name",
        target_index="gene_name",
        edge_attr=["miRNA Interaction Site", "Transcript ID"],
        directed=True,
        relabel_nodes=None,
        npartitions=0,
    ):
        if file_resources is None:
            file_resources = {}
            file_resources["miRNA_binding_sites.txt"] = os.path.join(
                path, "miRNA_binding_sites.txt")
            file_resources["general_information.txt"] = os.path.join(
                path, "general_information.txt")

        super(lncRNome, self).__init__(
            path,
            file_resources,
            source_col_name,
            target_col_name,
            source_index,
            target_index,
            edge_attr,
            directed,
            relabel_nodes,
            npartitions,
        )

    def load_network(
        self,
        file_resources,
        source_col_name,
        target_col_name,
        edge_attr,
        directed,
        filters,
    ):
        df = pd.read_table(self.file_resources["miRNA_binding_sites.txt"],
                           header=0)
        logging.info(f"{self.name()}, {df.columns.tolist()}")

        df["Binding miRNAs"] = df["Binding miRNAs"].str.lower()
        df["Binding miRNAs"] = df["Binding miRNAs"].str.replace(
            "-3p.*|-5p.*", "")

        lncRNome_miRNA_binding_sites_network = nx.from_pandas_edgelist(
            df,
            source=source_col_name,
            target=target_col_name,
            edge_attr=edge_attr,
            create_using=nx.DiGraph(),
        )

        return lncRNome_miRNA_binding_sites_network

    def load_dataframe(self, file_resources, npartitions=None):
        return pd.read_table(
            self.file_resources["general_information.txt"],
            header=0,
            usecols=[
                "Gene Name",
                "Transcript Name",
                "Transcript Type",
                "Location",
                "Strand",
            ],
        )


class NPInter(Interactions):
    def __init__(
        self,
        path="http://bigdata.ibp.ac.cn/npinter4/download/",
        file_resources=None,
        source_col_name="ncName",
        target_col_name="tarName",
        source_index="gene_name",
        target_index="gene_name",
        edge_attr=["tarType", "tissueOrCell", "tag", "level"],
        filters=None,
        directed=True,
        relabel_nodes=None,
        verbose=False,
    ):
        if file_resources is None:
            file_resources = {}
            file_resources["interaction_NPInterv4.expr.txt"] = os.path.join(
                path, "file/interaction_NPInterv4.expr.txt.gz")

        super(NPInter, self).__init__(
            path=path,
            file_resources=file_resources,
            source_col_name=source_col_name,
            target_col_name=target_col_name,
            source_index=source_index,
            target_index=target_index,
            edge_attr=edge_attr,
            filters=filters,
            directed=directed,
            relabel_nodes=relabel_nodes,
            verbose=verbose,
        )

    def load_network(
        self,
        file_resources,
        source_col_name,
        target_col_name,
        edge_attr,
        directed,
        filters,
    ):
        df = pd.read_table(file_resources["interaction_NPInterv4.expr.txt"],
                           header=0,
                           na_values=["-"])
        logging.info(f"{self.name()}, {df.columns.tolist()}")
        df["ncName"] = df["ncName"].str.upper()
        df["ncName"] = df["ncName"].str.strip("LNCRNA-")
        df["ncName"] = df["ncName"].str.replace("MALAT-1", "MALAT1")
        df["ncName"] = df["ncName"].str.replace("^MIR-",
                                                "hsa-mir-",
                                                regex=True)
        df["ncName"] = df["ncName"].str.replace("^MICRORNA-",
                                                "hsa-mir-",
                                                regex=True)

        df["tarName"] = df["tarName"].str.upper()

        df = self.filter_values(df, filters)

        lncRNome_miRNA_binding_sites_network = nx.from_pandas_edgelist(
            df,
            source=source_col_name,
            target=target_col_name,
            edge_attr=edge_attr,
            create_using=nx.DiGraph() if directed else nx.Graph(),
        )

        return lncRNome_miRNA_binding_sites_network


class StarBase(Interactions):
    def __init__(
        self,
        path,
        file_resources,
        source_col_name="geneName",
        target_col_name="pairGeneName",
        source_index="gene_name",
        target_index="gene_name",
        min_interactionNum=1,
        min_expNum=1,
        edge_attr=None,
        directed=True,
        relabel_nodes=None,
        npartitions=0,
    ):
        if file_resources is None:
            file_resources = {}
            file_resources[
                "starbase_3.0_lncrna_rna_interactions.csv"] = os.path.join(
                    path, "starbase_3.0_lncrna_rna_interactions.csv")
        self.min_interactionNum = min_interactionNum
        self.min_expNum = min_expNum
        super(StarBase, self).__init__(
            path,
            file_resources,
            source_col_name,
            target_col_name,
            source_index,
            target_index,
            edge_attr,
            directed,
            relabel_nodes,
            npartitions,
        )

    def load_network(
        self,
        file_resources,
        source_col_name,
        target_col_name,
        edge_attr,
        directed,
        filters,
    ):
        df = pd.read_csv(
            self.file_resources["starbase_3.0_lncrna_rna_interactions.csv"],
            header=0)

        df.loc[df["pairGeneType"] == "miRNA", "pairGeneName"] = df[
            df["pairGeneType"] == "miRNA"]["pairGeneName"].str.lower()
        df.loc[df["pairGeneType"] == "miRNA", "pairGeneName"] = df[
            df["pairGeneType"] == "miRNA"]["pairGeneName"].str.replace(
                "-3p.*|-5p.*", "")
        df = df[df["interactionNum"] >= self.min_interactionNum]
        df = df[df["expNum"] >= self.min_expNum]

        self.starBase_RNA_RNA_network = nx.from_pandas_edgelist(
            df,
            source=source_col_name,
            target=target_col_name,
            edge_attr=["interactionNum"],
            create_using=nx.DiGraph(),
        )
        return self.starBase_RNA_RNA_network


class MiRTarBase(Interactions):
    def __init__(
        self,
        path="http://mirtarbase.mbc.nctu.edu.tw/cache/download/7.0/",
        file_resources=None,
        source_col_name="miRNA",
        target_col_name="Target Gene",
        source_index="transcript_name",
        target_index="gene_name",
        edge_attr=None,
        filters={"Species (Target Gene)": "Homo sapiens"},
        directed=True,
        relabel_nodes=None,
        strip_mirna_name=False,
    ):
        if edge_attr is None:
            edge_attr = ["Support Type"]
        self.strip_mirna_name = strip_mirna_name

        if file_resources is None:
            file_resources = {}
            file_resources["miRTarBase_MTI.xlsx"] = os.path.join(
                path, "miRTarBase_MTI.xlsx")

        super(MiRTarBase, self).__init__(
            path=path,
            file_resources=file_resources,
            source_col_name=source_col_name,
            target_col_name=target_col_name,
            source_index=source_index,
            target_index=target_index,
            edge_attr=edge_attr,
            filters=filters,
            directed=directed,
            relabel_nodes=relabel_nodes,
        )

    def load_network(
        self,
        file_resources,
        source_col_name,
        target_col_name,
        edge_attr,
        directed,
        filters,
    ):
        df = pd.read_excel(self.file_resources["miRTarBase_MTI.xlsx"])
        logging.info(f"{self.name()}, {df.columns.tolist()}")

        df = self.filter_values(df, filters)

        if self.strip_mirna_name:
            df["miRNA"] = df["miRNA"].str.lower()
            df["miRNA"] = df["miRNA"].str.replace("-3p.*|-5p.*", "")

        mir_target_network = nx.from_pandas_edgelist(
            df,
            source=source_col_name,
            target=target_col_name,
            edge_attr=edge_attr,
            create_using=nx.DiGraph() if directed else nx.Graph(),
        )
        return mir_target_network


class TargetScan(Interactions, Dataset):
    def __init__(
        self,
        path,
        file_resources=None,
        source_col_name="MiRBase ID",
        target_col_name="Gene Symbol",
        source_index="transcript_name",
        target_index="transcript_name",
        edge_attr=["tissue", "positive_negative"],
        directed=True,
        relabel_nodes=None,
        species=9606,
        strip_mirna_name=False,
    ):
        if edge_attr is None:
            edge_attr = ["tissue", "positive_negative"]
        self.strip_mirna_name = strip_mirna_name
        self.species = species
        if file_resources is None:
            file_resources = {}
            file_resources["miR_Family_Info.txt"] = os.path.join(
                path, "miR_Family_Info.txt")
            file_resources[
                "Predicted_Targets_Info.default_predictions.txt"] = os.path.join(
                    path, "Predicted_Targets_Info.default_predictions.txt")

        super(TargetScan, self).__init__(
            path=path,
            file_resources=file_resources,
            source_col_name=source_col_name,
            target_col_name=target_col_name,
            source_index=source_index,
            target_index=target_index,
            edge_attr=edge_attr,
            directed=directed,
            relabel_nodes=relabel_nodes,
        )

    def load_network(
        self,
        file_resources,
        source_col_name,
        target_col_name,
        edge_attr,
        directed,
        filters,
    ):
        self.df = self.process_miR_family_info_table(file_resources,
                                                     self.species)
        interactions_df = self.process_interactions_table(
            file_resources, self.df, self.species)
        logging.info(f"{self.name()}, {interactions_df.columns.tolist()}")

        mir_target_network = nx.from_pandas_edgelist(
            interactions_df,
            source=source_col_name,
            target=target_col_name,
            edge_attr=edge_attr,
            create_using=nx.DiGraph() if directed else nx.Graph(),
        )
        return mir_target_network

    def process_miR_family_info_table(self, file_resources, species=None):
        miR_Family_Info_df = pd.read_table(
            file_resources["miR_Family_Info.txt"], delimiter="\t")

        if species:
            miR_Family_Info_df = miR_Family_Info_df[
                miR_Family_Info_df["Species ID"] == species]

        # Standardize MiRBase ID to miRNA names obtained from RNA-seq hg19
        if self.strip_mirna_name:
            miR_Family_Info_df["MiRBase ID"] = miR_Family_Info_df[
                "MiRBase ID"].str.lower()
            miR_Family_Info_df["MiRBase ID"] = miR_Family_Info_df[
                "MiRBase ID"].str.replace("-3p.*|-5p.*", "")

        miR_Family_Info_df.drop_duplicates(inplace=True)
        miR_Family_Info_df = miR_Family_Info_df.filter(
            items=[
                "miR family",
                "MiRBase ID",
                "Seed+m8",
                "Mature sequence",
                "Family Conservation?",
                "MiRBase Accession",
            ],
            axis="columns",
        )
        miR_Family_Info_df["MiRBase ID"] = miR_Family_Info_df[
            "MiRBase ID"].astype(str)
        return miR_Family_Info_df

    def process_interactions_table(self, file_resources, family_to_miR_df,
                                   species):
        """
        This functions joins the interactions data table between miR Family and targets, and
        Args:
            file_resources:
            family_to_miR_df:
            species:

        Returns:

        """
        # Load data frame from file
        family_interactions_df = pd.read_table(
            file_resources["Predicted_Targets_Info.default_predictions.txt"],
            delimiter="\t",
            low_memory=True,
        )

        # Select only homo sapiens miRNA-target pairs
        if species:
            family_interactions_df = family_interactions_df[
                family_interactions_df["Species ID"] == species]

        family_interactions_df = family_interactions_df.filter(
            items=["miR Family", "Gene Symbol"], axis="columns")
        family_to_miR_df = family_to_miR_df.filter(
            items=["miR family", "MiRBase ID"], axis="columns")
        family_to_miR_df.rename(columns={"miR family": "miR Family"},
                                inplace=True)

        # map miRBase ID names to miR Family
        # family_interactions_df = pd.merge(family_interactions_df, family_to_miR_df, how='outer', on="miR Family")

        family_to_miR_df.set_genes_index("miR Family", inplace=True)
        family_interactions_df.set_genes_index("miR Family", inplace=True)
        mir_interactions_df = family_interactions_df.join(
            family_to_miR_df, how="outer", on="miR Family").reset_index()

        # Standardize MiRBase ID to miRNA names obtained from RNA-seq hg19
        if self.strip_mirna_name:
            mir_interactions_df["MiRBase ID"] = mir_interactions_df[
                "MiRBase ID"].str.lower()
            mir_interactions_df["MiRBase ID"] = mir_interactions_df[
                "MiRBase ID"].str.replace("-3p.*|-5p.*", "")

        return mir_interactions_df
