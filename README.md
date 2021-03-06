<p align="center">
  <img src="https://github.com/BioMeCIS-Lab/OpenOmics/raw/master/openomics_web/assets/openomics_logo.png" max-height="200">
</p>

[![PyPI version](https://badge.fury.io/py/openomics.svg)](https://badge.fury.io/py/openomics)
[![Documentation Status](https://readthedocs.org/projects/openomics/badge/?version=latest)](https://openomics.readthedocs.io/en/latest/?badge=latest)
[![pyOpenSci](https://tinyurl.com/y22nb8up)](https://github.com/pyOpenSci/software-review/issues/31)
[![status](https://joss.theoj.org/papers/aca43e3c2989a803b514faef72dd3294/status.svg)](https://joss.theoj.org/papers/aca43e3c2989a803b514faef72dd3294)
[![DOI](https://zenodo.org/badge/125549505.svg)](https://zenodo.org/badge/latestdoi/125549505)
[![OpenOmics](https://github.com/BioMeCIS-Lab/OpenOmics/actions/workflows/python-package.yml/badge.svg?branch=master)](https://github.com/BioMeCIS-Lab/OpenOmics/actions/workflows/python-package.yml)
[![codecov](https://codecov.io/gh/BioMeCIS-Lab/OpenOmics/branch/master/graph/badge.svg?token=6N1UZ27MPH)](https://codecov.io/gh/BioMeCIS-Lab/OpenOmics)

This Python package provide a series of tools to integrate and explore the genomics, transcriptomics, proteomics, and
clinical data (aka multi-omics data). With interfaces to popular annotation databases and scalable data-frame manipulation tools, OpenOmics facilitates the common
data wrangling tasks when preparing data for RNA-seq bioinformatics analysis.

Documentation ([Latest](https://openomics.readthedocs.io/en/latest/)
| [Stable](https://openomics.readthedocs.io/en/stable/))
| [OpenOmics at a glance](https://openomics.readthedocs.io/en/latest/usage/getting-started.html)

## Features
OpenOmics assist in integration of heterogeneous multi-omics bioinformatics data. The library provides a Python API as well as an interactive Dash web interface.
It features support for:
- Genomics, Transcriptomics, Proteomics, and Clinical data.
- Harmonization with 20+ popular annotation, interaction, disease-association databases.

OpenOmics also has an efficient data pipeline that bridges the popular data manipulation Pandas library and Dask distributed processing to address the following use cases:

- Providing a standard pipeline for dataset indexing, table joining and querying, which are transparent and customizable
  for end-users.
- Providing Efficient disk storage for large multi-omics dataset with Parquet data structures.
- Integrating various data types including interactions and sequence data, then exporting to NetworkX graphs or data generators for down-stream machine learning.
- Accessible by both developers and scientists with a Python API that works seamlessly with an external Galaxy tool interface or the built-in Dash web interface (WIP).


## Installation via pip:

```
$ pip install openomics
```

## 

## Citations
The journal paper for this scientific package is currently being reviewed. In the meanwhile, the current package version can be cited with:

    # BibTeX
    @software{nhat_tran_2021_4731011,
      author       = {Nhat Tran and
                      Jean Gao},
      title        = {{OpenOmics: A bioinformatics API and web-app to 
                       integrate multi-omics datasets and interface with
                       public databases.}},
      month        = apr,
      year         = 2021,
      publisher    = {Zenodo},
      version      = {v0.8.8},
      doi          = {10.5281/zenodo.4731011},
      url          = {https://doi.org/10.5281/zenodo.4731011}
    }


## Credits

This package was created with [Cookiecutter](https://github.com/audreyr/cookiecutter) and the [pyOpenSci/cookiecutter-pyopensci](https://github.com/pyOpenSci/cookiecutter-pyopensci) project template, based off [audreyr/cookiecutter-pypackage](https://github.com/audreyr/cookiecutter-pypackage).
