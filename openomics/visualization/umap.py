import numpy as np
import plotly.express as px
import umap


def d3_umap(X, y_km, heat=None):
    """
    Args:
        X:
        y_km:
        heat:
    """
    reducer = umap.UMAP(random_state=1234, n_components=3)
    X_embedded = reducer.fit_transform(X)
    node_colors = get_node_colormap(y_km)
    x, y, z = X_embedded[:, 0], X_embedded[:, 1], X_embedded[:, 2]

    fig = px.scatter_3d(x=x, y=y, z=z, color=node_colors)
    fig.show()
    return reducer


def get_node_colormap(node_label):
    """
    Args:
        node_label:
    """
    if type(node_label) == list:
        node_labels = node_label
        sorted_node_labels = sorted(set(node_labels), reverse=True)
        colors = np.linspace(0, 1, len(sorted_node_labels))
        node_colormap = {
            f: colors[sorted_node_labels.index(f)]
            for f in set(node_labels)
        }
        node_colors = [
            node_colormap[n] if n in node_colormap.keys() else None
            for n in node_labels
        ]

    elif node_label.dtype == "object":
        node_labels = node_label.str.split("|", expand=True)[0]
        sorted_node_labels = sorted(node_labels.unique(), reverse=True)
        colors = np.linspace(0, 1, len(sorted_node_labels))
        node_colormap = {
            f: colors[sorted_node_labels.index(f)]
            for f in node_labels.unique()
        }
        node_colors = [
            node_colormap[n] if n in node_colormap.keys() else None
            for n in node_labels
        ]

    elif node_label.dtype == "float":
        node_labels = node_label.values
        node_colormap = None
        node_colors = [n / node_labels.max() for n in node_labels]
    return node_colors
