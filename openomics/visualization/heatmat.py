import pandas as pd
import plotly.graph_objects as go


def heatmap(table,
            file_output=None,
            title=None,
            autosize=True,
            width=800,
            height=1000):
    """
    Args:
        table:
        file_output:
        title:
        autosize:
        width:
        height:
    """
    if type(table.columns) == pd.MultiIndex:
        columns = table.columns.to_series().apply(
            lambda x: "{0}-{1}".format(*x))
    else:
        columns = table.columns
    fig = go.Figure(data=go.Heatmap(
        z=table,
        x=columns,
        y=table.index,
        hoverongaps=False,
    ))
    fig.update_layout(
        title=title,
        autosize=autosize,
        width=width,
        height=height,
    )
    if file_output:
        fig.write_image(file_output)

    return fig
