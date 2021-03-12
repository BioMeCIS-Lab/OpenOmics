import numpy as np
import pandas as pd


def concat_uniques(series: pd.Series):
    """An aggregation custom function to be applied to each column of a groupby
    Args:
        series (pd.Series):
    """
    series_str = series.dropna().astype(str)
    if not series_str.empty:
        return "|".join(series_str.unique())
    else:
        return None


def concat(series: pd.Series):
    """
    Args:
        series (pd.Series):
    """
    series = series.dropna().astype(str)
    if not series.empty:
        return "|".join(series)
    else:
        return None


def drop_duplicate_columns(df):
    """
    Args:
        df:
    """
    _, i = np.unique(df.columns, return_index=True)
    df = df.iloc[:, i]
    return df


def slice_adj(adj, node_list: list, nodes_A, nodes_B=None):
    """
    Args:
        adj:
        node_list (list):
        nodes_A:
        nodes_B:
    """
    if nodes_B is None:
        idx = [node_list.index(node) for node in nodes_A]
        return adj[idx, :][:, idx]
    else:
        idx_A = [node_list.index(node) for node in nodes_A]
        idx_B = [node_list.index(node) for node in nodes_B]
        return adj[idx_A, :][:, idx_B]
