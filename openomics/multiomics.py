from typing import List, Dict, Union

from openomics import backend as pd
from .clinical import (
    ClinicalData,
    HISTOLOGIC_SUBTYPE,
    PATHOLOGIC_STAGE,
    TUMOR_NORMAL,
    PREDICTED_SUBTYPE,
)
from .genomics import SomaticMutation, CopyNumberVariation, DNAMethylation
from .imageomics import WholeSlideImage
from .proteomics import Protein
from .transcriptomics import MessengerRNA, MicroRNA, LncRNA, ExpressionData


class MultiOmics:
    def __init__(self, cohort_name):
        """Load all multi-omics data from a given cohort_folder path.

        Args:
            cohort_name (str): the clinical cohort name
        """
        self._cohort_name = cohort_name
        self._omics = []

        # This is a data dictionary accessor to retrieve individual -omic data
        self.data = {}

    def add_omic(self,
                 omic_data: ExpressionData,
                 initialize_annotations: bool = True):
        """Adds an omic object to the Multiomics such that the samples in omic
        matches the samples existing in the other omics.

        Args:
            omic_data (ExpressionData): The omic to add, e.g., MessengerRNA,
                MicroRNA, LncRNA, etc.
            initialize_annotations (bool): default True. If true, initializes
                the annotation dataframe in the omic object
        """
        self.__dict__[omic_data.name()] = omic_data

        if omic_data.name not in self._omics:
            self._omics.append(omic_data.name())

        # dictionary as data accessor to the expression data
        self.data[omic_data.name()] = omic_data.expressions

        # Initialize annotation
        if initialize_annotations:
            omic_data.initialize_annotations(gene_list=None,
                                             index=omic_data.gene_index)

        print(
            omic_data.name(),
            self.data[omic_data.name()].shape if hasattr(
                self.data[omic_data.name()], "shape") else ": None",
            ", indexed by:",
            omic_data.annotations.index.name,
        )

    def add_clinical_data(self, clinical_data=None):
        """
        Args:
            clinical_data:
        """
        if type(clinical_data) == ClinicalData:
            self.clinical = clinical_data
        else:
            self.clinical = ClinicalData(clinical_data)

        self.data["PATIENTS"] = self.clinical.patient
        if hasattr(self.clinical, "biospecimen"):
            self.data["BIOSPECIMENS"] = self.clinical.biospecimen
        if hasattr(self.clinical, "drugs"):
            self.data["DRUGS"] = self.clinical.drugs

    def get_omics_list(self):
        return self._omics

    def __getitem__(self, item):
        # type: (str) -> object
        """This function allows the MultiOmicData class objects to access
        individual omics by a dictionary lookup, e.g. openomics["MicroRNA"]

        Args:
            item (str): a string of the class name
        """
        if item.lower() == MessengerRNA.name().lower():
            return self.__getattribute__(MessengerRNA.name())

        elif item.lower() == MicroRNA.name().lower():
            return self.__getattribute__(MicroRNA.name())

        elif item.lower() == LncRNA.name().lower():
            return self.__getattribute__(LncRNA.name())

        elif item.lower() == WholeSlideImage.name().lower():
            return self.__getattribute__(WholeSlideImage.name())

        elif item.lower() == SomaticMutation.name().lower():
            return self.__getattribute__(SomaticMutation.name())

        elif item.lower() == CopyNumberVariation.name().lower():
            return self.__getattribute__(CopyNumberVariation.name())

        elif item.lower() == DNAMethylation.name().lower():
            return self.__getattribute__(DNAMethylation.name())

        elif item.lower() == Protein.name().lower():
            return self.__getattribute__(Protein.name())

        elif item.lower() == "patients":
            return self.clinical.patient
        elif item.lower() == "samples":
            return self.clinical.samples
        elif item.lower() == "drugs":
            return self.clinical.drugs
        else:
            raise Exception(
                'String accessor must be one of {"MessengerRNA", "MicroRNA", "LncRNA", "Protein", etc.}'
            )

    def remove_duplicate_genes(self):
        """Removes duplicate genes between any omics such that the gene index
        across all omics has no duplicates.
        """
        for omic_A in self._omics:
            for omic_B in self._omics:
                if omic_A != omic_B:
                    self.__getattribute__(omic_A).drop_genes(
                        set(self.__getattribute__(omic_A).get_genes_list())
                        & set(self.__getattribute__(omic_B).get_genes_list()))

    def build_samples(self, agg_by="union"):
        """Running this function will build a dataframe for all samples across
        the different omics (either by a union or intersection). Then,

        Args:
            agg_by (str): ["union", "intersection"]
        """
        if len(self._omics
               ) < 1:  # make sure at least one ExpressionData present
            print(
                "build_samples() does nothing. Must add at least one omic to this MultiOmics object."
            )

        all_samples = pd.Index([])
        for omic in self._omics:
            if agg_by == "union":
                all_samples = all_samples.union(self.data[omic].index)
            elif agg_by == "intersection":
                all_samples = all_samples.intersection(self.data[omic].index)

        if hasattr(self, "clinical"):
            self.clinical.build_clinical_samples(all_samples)
            self.data["SAMPLES"] = self.clinical.samples.index
        else:
            self.data["SAMPLES"] = all_samples

    def __dir__(self):
        return list(self.data.keys())

    def match_samples(self, omics) -> pd.Index:
        """Return the index of bcr_sample_barcodes of the intersection of
        samples from all modalities

        Args:
            omics: An array of modalities

        Returns:
            matched_sapmles: An pandas Index list
        """
        # TODO check that for single modalities, this fetch all patients
        matched_samples = self.data[omics[0]].index.copy()

        for omic in omics:
            matched_samples = matched_samples.join(self.data[omic].index,
                                                   how="inner")

        return matched_samples

    def load_data(
        self,
        omics,
        target=["pathologic_stage"],
        pathologic_stages=None,
        histological_subtypes=None,
        predicted_subtypes=None,
        tumor_normal=None,
        samples_barcode=None,
    ):
        # type: (Union[List[str], str], List[str], List[str], List[str], List[str], List[str], List[str]) -> (Dict[str, pd.DataFrame], pd.DataFrame)
        """
        Args:
            omics (list): A list of the data modalities to load. Default "all"
                to select all modalities
            target (list): The clinical data fields to include in the
            pathologic_stages (list): Only fetch samples having certain stages
                in their corresponding patient's clinical data. For instance,
                ["Stage I", "Stage II"] will only fetch samples from Stage I and
                Stage II patients. Default is [] which fetches all pathologic
                stages.
            histological_subtypes: A list specifying the histological subtypes
                to fetch. Default is [] which fetches all histological sybtypes.
            predicted_subtypes: A list specifying the histological subtypes to
                fetch. Default is [] which fetches all histological sybtypes.
            tumor_normal: ["Tumor", "Normal"]. Default is [], which fetches all
                tumor or normal sample types.
            samples_barcode: A list of sample's barcode. If not None, only fetch
                data with matching samples provided in this list.

        Returns:
            (X, y): Returns X, a dictionary containing the multiomics data that
            have data
        """
        if omics == "all" or omics is None:
            omics = self._omics

        matched_samples = self.match_samples(omics)

        if samples_barcode is not None:
            matched_samples = samples_barcode

        if hasattr(self, "clinical") and isinstance(self.clinical,
                                                    ClinicalData):
            # Build targets clinical data
            y = self.get_patients_clinical(matched_samples)

            # Select only samples with certain cancer stage or subtype
            if pathologic_stages:
                y = y[y[PATHOLOGIC_STAGE].isin(pathologic_stages)]
            if histological_subtypes:
                y = y[y[HISTOLOGIC_SUBTYPE].isin(histological_subtypes)]
            if predicted_subtypes:
                y = y[y[PREDICTED_SUBTYPE].isin(predicted_subtypes)]
            if tumor_normal:
                y = y[y[TUMOR_NORMAL].isin(tumor_normal)]

            # Filter y target column labels
            y = y.filter(target)
            y.dropna(axis=0, inplace=True)
            matched_samples = y.index
        else:
            y = None

        # Build expression matrix for each omic, indexed by matched_samples
        X_multiomics = {}
        for omic in omics:
            X_multiomics[omic] = self.data[omic].loc[
                matched_samples, self[omic].get_genes_list()]

        return X_multiomics, y

    def get_patients_clinical(self, matched_samples):
        """Fetch patient's clinical data for each given samples barcodes in the
        matched_samples

        Returns
            samples_index: Index of samples

        Args:
            matched_samples: A list of sample barcodes
        """
        return self.data["SAMPLES"].reindex(matched_samples)

    def print_sample_sizes(self):
        for omic in self.data:
            print(
                omic,
                self.data[omic].shape
                if hasattr(self.data[omic], "shape") else "Didn't import data",
            )

    def add_subtypes_to_patients_clinical(self, dictionary):
        """This function adds a "predicted_subtype" field to the patients
        clinical data. For instance, patients were classified into subtypes
        based on their expression profile using k-means, then, to use this
        function, do:

        add_subtypes_to_patients_clinical(dict(zip(<list of patient
        barcodes>, <list of corresponding patient's subtypes>)))

        Adding a field to the patients clinical data allows openomics to
        query the patients data through the .load_data(predicted_subtypes=[])
        parameter,

        Args:
            dictionary: A dictionary mapping patient's barcode to a subtype
        """
        self.data["PATIENTS"] = self.data["PATIENTS"].assign(
            predicted_subtype=self.data["PATIENTS"][
                self.clinical.patient_column].map(dictionary))
