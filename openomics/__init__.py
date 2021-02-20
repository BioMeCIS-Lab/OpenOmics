from __future__ import print_function, division, absolute_import

import json
import logging
import os
import sys

import astropy
import dask.dataframe as dd
import pandas as pd

"""Top-level package for openomics."""

__author__ = """Nhat (Jonny) Tran"""
__email__ = 'nhat.tran@mavs.uta.edu'
__version__ = '0.8.5'


# Initialize configurations
this = sys.modules[__name__]
this.config = {}

# Set pandas backend
this.config["backend"] = pd

# Set cache download directory
this.config["cache_dir"] = astropy.config.get_cache_dir(this.__name__)
logging.info("Cache directory is", this.config["cache_dir"])

user_conf_path = "~/.openomics/conf.json"

# Initialize user configuration file at ~/.openomics/conf.json
if not os.path.exists(user_conf_path):
    if not os.path.exists("~/.openomics"):
        os.makedirs("~/.openomics")

    if not os.path.exists(user_conf_path):
        base_config = {}
        base_config['cache_dir'] = astropy.config.get_cache_dir(this.__name__)

        with open(user_conf_path, 'w', encoding='utf-8') as file:
            json.dump(base_config, fp=file, indent=4)

# Read configuration from ~/.openomics/conf.json
if os.path.exists(user_conf_path):
    with open(user_conf_path, 'a', encoding='utf-8') as file:
        user_config = json.load(fp=file)

    if user_config:
        for p in user_config['database']:
            this.config.update(p)

from . import database, utils

from .transcriptomics import (
    Expression, MessengerRNA, MicroRNA, LncRNA,
)

from .genomics import (
    SomaticMutation, DNAMethylation, CopyNumberVariation
)

from .proteomics import (
    Protein
)

from .clinical import ClinicalData

from .multiomics import (
    MultiOmics
)


def set_backend(new):
    """
    Args:
        new:
    """
    assert new in ["dask", "pandas"]

    if new == "dask":
        this.config["backend"] = dd
    else:
        this.config["backend"] = pd


def set_cache_dir(path, delete_temp=False):
    if not os.path.exists(path):
        raise NotADirectoryError(path)

    this.config["cache_dir"] = path
    astropy.config.set_temp_cache(path=path, delete=delete_temp)
    logging.info("Cache directory is", this.config["cache_dir"])
