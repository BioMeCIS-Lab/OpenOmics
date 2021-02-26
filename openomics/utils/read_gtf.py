# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

from __future__ import print_function, division, absolute_import

import logging
from collections import OrderedDict
from os.path import exists

import dask.dataframe as dd
import numpy as np
import pandas as pd
from six import string_types
from six.moves import intern


def expand_attribute_strings(attribute_strings,
                             quote_char='"',
                             nan_value=None,
                             usecols=None):
    """The last column of GTF has a variable number of key value pairs of the
    format: "key1 value1; key2 value2;" Parse these into a dictionary mapping
    each key onto a list of values, where the value is None for any row where
    the key was missing.

    Args:
        attribute_strings (list of str):
        quote_char (str): Quote character to remove from values
        nan_value (any): If an attribute is missing from a row, give it this
            value.
        usecols (list of str or None): If not None, then only expand columns
            included in this set, otherwise use all columns.
    """
    n = len(attribute_strings)

    extra_columns = {}
    column_order = []

    #
    # SOME NOTES ABOUT THE BIZARRE STRING INTERNING GOING ON BELOW
    #
    # While parsing millions of repeated strings (e.g. "gene_id" and "TP53"),
    # we can save a lot of memory by making sure there's only one string
    # object per unique string. The canonical way to do this is using
    # the 'intern' function. One problem is that Py2 won't let you intern
    # unicode objects, so to get around this we call intern(str(...)).
    #
    # It also turns out to be faster to check interned strings ourselves
    # using a local dictionary, hence the two dictionaries below
    # and pair of try/except blocks in the loop.
    column_interned_strings = {}
    value_interned_strings = {}

    for (i, attribute_string) in enumerate(attribute_strings):
        for kv in attribute_string.split(";"):
            # We're slicing the first two elements out of split() because
            # Ensembl release 79 added values like:
            #   transcript_support_level "1 (assigned to previous version 5)";
            # ...which gets mangled by splitting on spaces.
            parts = kv.strip().split(" ", 2)[:2]

            if len(parts) != 2:
                continue

            column_name, value = parts

            try:
                column_name = column_interned_strings[column_name]
            except KeyError:
                column_name = intern(str(column_name))
                column_interned_strings[column_name] = column_name

            if usecols is not None and column_name not in usecols:
                continue

            try:
                column = extra_columns[column_name]
            except KeyError:
                column = [nan_value] * n
                extra_columns[column_name] = column
                column_order.append(column_name)

            value = (value.replace(quote_char, "")
                     if value.startswith(quote_char) else value)

            try:
                value = value_interned_strings[value]
            except KeyError:
                value = intern(str(value))
                value_interned_strings[value] = value

            # if an attribute is used repeatedly then
            # keep track of all its values in a list
            old_value = column[i]
            if old_value is nan_value:
                column[i] = value
            else:
                column[i] = "%s,%s" % (old_value, value)

    return OrderedDict((column_name, extra_columns[column_name])
                       for column_name in column_order)


REQUIRED_COLUMNS = [
    "seqname",
    "source",
    "feature",
    "start",
    "end",
    "score",
    "strand",
    "frame",
    "attribute",
]

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def parse_gtf(
    filepath_or_buffer,
    chunksize=1024 * 1024,
    features=None,
    intern_columns=["seqname", "source", "strand", "frame"],
    fix_quotes_columns=["attribute"],
):
    """
    Args:
        filepath_or_buffer (str or buffer object):
        chunksize (int): Default 1048576.
        features (set or None): Drop entries which aren't one of these features
        intern_columns (list): These columns are short strings which should be
            interned
        fix_quotes_columns (list): Most commonly the 'attribute' column which
            had broken quotes on some Ensembl release GTF files.
    """
    if features is not None:
        features = set(features)

    def parse_frame(s):
        if s == ".":
            return 0
        else:
            return int(s)

    # GTF columns:
    # 1) seqname: str ("1", "X", "chrX", etc...)
    # 2) source : str
    #      Different versions of GTF use second column as of:
    #      (a) gene biotype
    #      (b) transcript biotype
    #      (c) the annotation source
    #      See: https://www.biostars.org/p/120306/#120321
    # 3) feature : str ("gene", "transcript", &c)
    # 4) start : int
    # 5) end : int
    # 6) score : float or "."
    # 7) strand : "+", "-", or "."
    # 8) frame : 0, 1, 2 or "."
    # 9) attribute : key-value pairs separated by semicolons
    # (see more complete description in docstring at top of file)

    chunk_iterator = pd.read_csv(
        filepath_or_buffer,
        sep="\t",
        comment="#",
        names=REQUIRED_COLUMNS,
        skipinitialspace=True,
        skip_blank_lines=True,
        error_bad_lines=True,
        warn_bad_lines=True,
        chunksize=chunksize,
        engine="c",
        dtype={
            "start": np.int64,
            "end": np.int64,
            "score": np.float32,
            "seqname": str,
        },
        na_values=".",
        converters={"frame": parse_frame},
    )

    dataframes = []
    try:
        for df in chunk_iterator:
            for intern_column in intern_columns:
                df[intern_column] = [intern(str(s)) for s in df[intern_column]]

            # compare feature strings after interning
            if features is not None:
                df = df[df["feature"].isin(features)]

            for fix_quotes_column in fix_quotes_columns:
                # Catch mistaken semicolons by replacing "xyz;" with "xyz"
                # Required to do this since the Ensembl GTF for Ensembl
                # release 78 has mistakes such as:
                #   gene_name = "PRAMEF6;" transcript_name = "PRAMEF6;-201"
                df[fix_quotes_column] = [
                    s.replace(';"', '"').replace(";-", "-")
                    for s in df[fix_quotes_column]
                ]
            dataframes.append(df)
    except Exception as e:
        raise e
        # raise Exception("ParsingError:" + str(e))

    df = pd.concat(dataframes)

    return df


def parse_gtf_dask(filepath_or_buffer,
                   npartitions=None,
                   compression=None,
                   features=None):
    """
    Args:
        filepath_or_buffer (str or buffer object):
        npartitions (int): Number of partitions for the dask dataframe. Default None.
        compression (str): Compression type to be passed into dask.dataframe.read_table(). Default None.
        features (set or None): Drop entries which aren't one of these features
    """
    if features is not None:
        features = set(features)

    def parse_frame(s):
        if s == ".":
            return 0
        else:
            return int(s)

    # GTF columns:
    # 1) seqname: str ("1", "X", "chrX", etc...)
    # 2) source : str
    #      Different versions of GTF use second column as of:
    #      (a) gene biotype
    #      (b) transcript biotype
    #      (c) the annotation source
    #      See: https://www.biostars.org/p/120306/#120321
    # 3) feature : str ("gene", "transcript", &c)
    # 4) start : int
    # 5) end : int
    # 6) score : float or "."
    # 7) strand : "+", "-", or "."
    # 8) frame : 0, 1, 2 or "."
    # 9) attribute : key-value pairs separated by semicolons
    # (see more complete description in docstring at top of file)

    # Uses Dask
    logging.debug("dask.datafame.read_table, file={}, compression={}".format(
        filepath_or_buffer, compression))
    dataframe = dd.read_table(
        filepath_or_buffer,
        sep="\t",
        compression=compression,
        blocksize=None,
        comment="#",
        names=REQUIRED_COLUMNS,
        skipinitialspace=True,
        skip_blank_lines=True,
        error_bad_lines=True,
        warn_bad_lines=True,
        # chunksize=chunksize,
        engine="c",
        dtype={
            "start": np.int64,
            "end": np.int64,
            "score": np.float32,
            "seqname": str,
        },
        na_values=".",
        converters={"frame": parse_frame},
    )

    return dataframe


def parse_gtf_and_expand_attributes(
    filepath_or_buffer,
    npartitions=None,
    compression=None,
    chunksize=1024 * 1024,
    restrict_attribute_columns=None,
    features=None,
):
    """Parse lines into column->values dictionary and then expand the
    'attribute' column into multiple columns. This expansion happens by
    replacing strings of semi-colon separated key-value values in the
    'attribute' column with one column per distinct key, with a list of values
    for each row (using None for rows where key didn't occur).

    Args:
        compression:
        filepath_or_buffer (str or buffer object):
        npartitions:
        chunksize (int):
        restrict_attribute_columns (list/set of str or None): If given, then only uses attribute columns.
        features (set or None): Ignore entries which don't correspond to one of the supplied features
    """
    if npartitions:
        df = parse_gtf_dask(
            filepath_or_buffer,
            npartitions=npartitions,
            compression=compression,
            features=features,
        )
        df = df.reset_index(drop=False)
        df = df.set_index("index")

        attribute_values = df.pop("attribute")

        for column_name, values in expand_attribute_strings(
                attribute_values, usecols=restrict_attribute_columns).items():
            series = dd.from_array(np.array(values, dtype=np.str))
            df[column_name] = series
    else:
        df = parse_gtf(filepath_or_buffer,
                       chunksize=chunksize,
                       features=features)

        attribute_values = df.pop("attribute")

        for column_name, values in expand_attribute_strings(
                attribute_values, usecols=restrict_attribute_columns).items():
            df[column_name] = values

    return df


def read_gtf(
    filepath_or_buffer,
    npartitions=None,
    compression=None,
    expand_attribute_column=True,
    infer_biotype_column=False,
    column_converters={},
    usecols=None,
    features=None,
    chunksize=1024 * 1024,
):
    """Parse a GTF into a dictionary mapping column names to sequences of
    values.

    Args:
        filepath_or_buffer (str or buffer object): Path to GTF file (may be gzip
            compressed) or buffer object such as StringIO
        npartitions (int): Number of partitions for the dask dataframe. Default None.
        compression (str): Compression type to be passed into dask.dataframe.read_table(). Default None.
        expand_attribute_column (bool): Replace strings of semi-colon separated
            key-value values in the 'attribute' column with one column per
            distinct key, with a list of values for each row (using None for
            rows where key didn't occur).
        infer_biotype_column (bool): Due to the annoying ambiguity of the second
            GTF column across multiple Ensembl releases, figure out if an older
            GTF's source column is actually the gene_biotype or
            transcript_biotype.
        column_converters (dict, optional): Dictionary mapping column names to
            conversion functions. Will replace empty strings with None and
            otherwise passes them to given conversion function.
        usecols (list of str or None): Restrict which columns are loaded to the
            give set. If None, then load all columns.
        features (set of str or None): Drop rows which aren't one of the
            features in the supplied set
        chunksize (int):
    """
    if isinstance(filepath_or_buffer,
                  string_types) and not exists(filepath_or_buffer):
        raise ValueError("GTF file does not exist: %s" % filepath_or_buffer)

    if expand_attribute_column:
        result_df = parse_gtf_and_expand_attributes(
            filepath_or_buffer,
            npartitions=npartitions,
            compression=compression,
            chunksize=chunksize,
            restrict_attribute_columns=usecols,
        )
    else:
        if npartitions:
            result_df = parse_gtf(filepath_or_buffer, features=features)
        else:
            result_df = parse_gtf_dask(
                filepath_or_buffer,
                npartitions=npartitions,
                features=features,
                compression=compression,
            )

    for column_name, column_type in list(column_converters.items()):
        result_df[column_name] = [
            column_type(string_value) if len(string_value) > 0 else None
            for string_value in result_df[column_name]
        ]

    # Hackishly infer whether the values in the 'source' column of this GTF
    # are actually representing a biotype by checking for the most common
    # gene_biotype and transcript_biotype value 'protein_coding'
    if infer_biotype_column:
        unique_source_values = set(result_df["source"])
        if "protein_coding" in unique_source_values:
            column_names = set(result_df.columns)
            # Disambiguate between the two biotypes by checking if
            # gene_biotype is already present in another column. If it is,
            # the 2nd column is the transcript_biotype (otherwise, it's the
            # gene_biotype)
            if "gene_biotype" not in column_names:
                logging.info(
                    "Using column 'source' to replace missing 'gene_biotype'")
                result_df["gene_biotype"] = result_df["source"]
            if "transcript_biotype" not in column_names:
                logging.info(
                    "Using column 'source' to replace missing 'transcript_biotype'"
                )
                result_df["transcript_biotype"] = result_df["source"]

    if usecols is not None:
        column_names = set(result_df.columns)
        valid_columns = [c for c in usecols if c in column_names]
        result_df = result_df[valid_columns]

    return result_df
