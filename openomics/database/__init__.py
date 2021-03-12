from __future__ import print_function, division, absolute_import

from .annotation import EnsemblGenes, RNAcentral, BioMartManager, GTEx, ProteinAtlas
from .base import Annotatable
from .disease import DisGeNet, HMDD, LncRNADisease, MalaCards
from .interaction import (
    LncBase,
    MiRTarBase,
    TargetScan,
    LncRNA2Target,
    BioGRID,
    GeneMania,
    lncRInter,
    lncRNome,
    STRING,
    NPInter,
)
from .ontology import GeneOntology
from .sequence import GENCODE, MirBase
