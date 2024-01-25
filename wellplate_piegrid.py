import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import string


def wellplate_piegrid(
    df: pd.DataFrame,
    data_columns,
    well_column,
    wells=96,
    colors=None,
    figsize=10,
    fontsize=14,
    legend_location: str = None,
    title: str = None,
    outfile=None,
    show=True,
):
    """
    Creates a grid of pie charts for well plates
    :param df: dataframe
    :param data_columns: all columns that make up the pie
    :param well_column: name of the well column that contains entries like A1, A2, ...
    :param wells: number of wells - 96 or 384 well plate
    :param colors: optional list of colors musst match up with number of data_columns
    :param figsize: figure size
    :param fontsize: font size for labels headers
    :param legend_location: None for no legend otherwise "right" or similar
    :param title:
    :param outfile: save to file
    :param show: show plot via plt.show()
    :return: the figure and axes objects
    """
    if colors is None:
        # extract colors from color map if not specified
        cmap = plt.get_cmap("tab10")
        colors = [cmap(i) for i in range(len(data_columns))]

    # Define headers for rows and columns
    # 96 well
    if wells == 96:
        cols = 12
        rows = 8
    else:  # 384 well
        cols = 24
        rows = 16

    # list of letters and numbers
    row_headers = list(string.ascii_uppercase)[:rows]
    col_headers = list(range(1, cols + 1))

    # Create a grid of pie charts
    fig, axs = plt.subplots(
        rows + 1, cols + 1, figsize=(figsize, figsize * cols / rows / 2.5)
    )
    # remove all placeholder plots
    for y in range(len(axs)):
        for x in range(len(axs[0])):
            axs[y, x].axis("off")

    # set headers
    for x in range(cols):
        axs[0, x + 1].text(
            0.5, 0.5, col_headers[x], ha="center", va="center", fontsize=fontsize
        )
    for y in range(rows):
        axs[y + 1, 0].text(
            0.5, 0.5, row_headers[y], ha="center", va="center", fontsize=fontsize
        )

    last_pie = None
    # create pies
    for i, (idx, row) in enumerate(df.iterrows()):
        well = row[well_column]
        y = ord(well[:1]) - 64  # A is 65
        x = int(well[1:])
        data = [row[group] for group in data_columns]
        axs[y, x].pie(data, labels=None, colors=colors)
        axs[y, x].axis("on")
        last_pie = axs[y, x]

    # legend
    if legend_location is not None:
        items = []
        for i, column in enumerate(data_columns):
            items.append(mpatches.Patch(color=colors[i], label=column))
        fig.legend(handles=items, loc=legend_location)

    if title is not None:
        plt.suptitle(title)

    if outfile is not None:
        plt.savefig(outfile, dpi=300)

    if show:
        plt.show()

    return fig, axs
