## Adapted from https://github.com/boecker-lab/myopic-mces-data/blob/main/umap_vis.py


import numpy as np
import pandas as pd
import plotly.express as px
from rdkit.Chem import PandasTools
from dash import dcc, html, Input, Output, no_update, Dash, callback
from rdkit.Chem import MolFromSmiles, Draw
import os
import pandas_utils as pu




COLOR_MAPPING = {'mcebio': '#3391C1',
                 'fdadrug': "#2F5697",
                 'nihnp': "#F08C2B",
                 'mcescaf': "#BF2C84",
                 'enam': "#06845D",
                 'enammol': "#FFE600",
                 'otavapep': "#4D2B6D",
                 "Missing": "#DADADA"
                 }

# for reverse proxy, add prefix to dash app settings
dash_settings = {}
try:
    for env in ["PROXY_PREFIX_ROUTES", "PROXY_PREFIX_REQUESTS", "PROXY_PREFIX_URL"]:
        if os.environ.get(env, "") != "":
            dash_settings.update(
                {
                    {
                        "PROXY_PREFIX_ROUTES": "routes_pathname_prefix",
                        "PROXY_PREFIX_REQUESTS": "requests_pathname_prefix",
                        "PROXY_PREFIX_URL": "url_base_pathname",
                    }[env]: os.environ[env]
                }
            )
except Exception as e:
    print("failed setting proxy setting in Dash", e)

app = Dash(__name__, title="TMAP_msnlib_chemical_space", **dash_settings)
server = app.server

umap_df = pu.read_dataframe(
    r"C:\git\msn_tree_library\data\acquisition_results\20240527_public_spectral_libraries_and_new_nist23_smiles_msnlib_divided_tmap_coord.tsv"
)
graph = dcc.Graph(id="umap-plot", clear_on_unhover=True, config=dict(scrollZoom=True))
libraries = [
    "msnlib",
    "massbankeu",
    "mona",
    "gnps",
    "mzcloud",
    "nist23",
    "open",
    "commercial",
    "available",
    "nihnp",
    "fdadrug",
    "mcescaf",
    "mcebio",
    "enammol",
    "enam",
    "otavapep"

]
dropdown = dcc.Dropdown(
    libraries,
    "msnlib",
    id="set-select",
    clearable=False,
)
tooltip = dcc.Tooltip(id="graph-tooltip")

app.layout = html.Div([dropdown, graph, tooltip], style=dict(width="200px"))


@app.callback(Output(graph, "figure"), Input(dropdown, "value"))
def update_figure(dataset):
    sorted = umap_df.sort_values(by=dataset)
    if(dataset == "msnlib"):

        fig = px.scatter(
            sorted,
            x="x",
            y="y",
            color='library',
            custom_data=["canonical_smiles", "entries", "availability"],
            template="simple_white",
            color_discrete_map=COLOR_MAPPING,
        )
        fig.update_traces(hoverinfo="none", hovertemplate=None)
        fig.update_traces(marker={"size": 3})
        fig.update_layout(
            dragmode="pan",
            width=970,
            height=800,
            legend={"itemsizing": "constant"},
            legend_title_text=None,
        )
    else:
        fig = px.scatter(
            sorted,
            x="x",
            y="y",
            color=dataset,
            custom_data=["canonical_smiles", "entries", "availability"],
            template="simple_white",
            color_discrete_map={False: "#aaaaaa", True: "#BF2C84"},
        )
        fig.update_traces(hoverinfo="none", hovertemplate=None)
        fig.update_traces(marker={"size": 3})
        fig.update_layout(
            dragmode="pan",
            width=970,
            height=800,
            legend={"itemsizing": "constant"},
            legend_title_text=None,
        )

    return fig


@app.callback(
    Output("graph-tooltip", "show"),
    Output("graph-tooltip", "bbox"),
    Output("graph-tooltip", "children"),
    Output("umap-plot", "extendData"),
    Input("umap-plot", "hoverData"),
)
def display_hover(hoverData):
    if hoverData is None:
        return False, no_update, no_update, no_update
    pt = hoverData["points"][0]
    if "customdata" not in pt:  # hovering over highlight
        return False, no_update, no_update, no_update
    bbox = pt["bbox"]
    # extract custom data
    smiles, entries, availability = pt["customdata"]
    ext_data = no_update
    data = Draw._moltoimg(MolFromSmiles(smiles), (200, 200), [], "", returnPNG=True)
    children = [
        html.Div(
            [
                html.H5(smiles),
                html.Img(
                    src=f"data:image/png;base64,{PandasTools._get_image(data)}",
                    alt="SMILES: " + smiles,
                    style=dict(
                        float="left",
                        margin="0px 15px 15px 0px",
                        height=200,
                        width=200,
                        border=2,
                    ),
                ),
                html.Div(
                    [
                        html.H5("Availability"),
                        html.Ul(
                            [
                                html.Li(c, style={"font-size": "10pt"})
                                for c in availability.split("_")
                            ]
                        ),
                    ]
                )
                # html.P(f'SMILES: {smiles}', style={'font-size': '6px'}),
            ],
            style={"width": "200px", "white-space": "normal"},
        )
    ]
    return True, bbox, children, ext_data


if __name__ == "__main__":
    for lib in libraries:
        fig = update_figure(lib)
        fig.update_layout(title=f"<b>{lib}</b>", title_x=0.5)
        # Define the directory and file name
        directory = r"C:\git\msn_tree_library\figures\public_msn_comparison"
        html = "{}.html".format(lib)
        pdf = "{}.pdf".format(lib)
        png = "{}.png".format(lib)
        # Ensure the directory exists
        if not os.path.exists(directory):
            os.makedirs(directory)
        # Define the full path
        html_file_path = os.path.join(directory, html)
        pdf_file_path = os.path.join(directory, pdf)
        png_file_path = os.path.join(directory, png)

        fig.write_html(html_file_path)
        fig.write_image(pdf_file_path)
        fig.write_image(png_file_path)
    app.run()
    app.write_html(r"C:\git\msn_library\figures\public_msn_comparison\all_libraries.html")
