import logging
import urllib.parse

import pandas as pd

from rest_utils import get_json_response, json_col

# NPAtlas:
# ORIGIN ORGANISM TYPE
NP_ATLAS_URL = "https://www.npatlas.org/api/v1/compounds/basicSearch?method=full&{}&threshold=0&orderby=npaid&ascending=true&limit=10"
NP_ATLAS_STRUCTURE_SEARCH_URL = "https://www.npatlas.org/api/v1/compounds/structureSearch?structure={}" \
                                "&type={}&method=sim&threshold=0.9999&skip=0&limit=10&stereo=false"
NP_ATLAS_SUFFIX = "_npatlas"


def np_atlas_url_by_inchikey(inchi_key):
    return NP_ATLAS_URL.format("inchikey=" + urllib.parse.quote(inchi_key))


def np_atlas_url_by_smiles(smiles):
    return NP_ATLAS_URL.format("smiles=" + urllib.parse.quote(smiles))


def np_atlas_url_similar_structure(structure, type="inchikey"):
    return NP_ATLAS_URL.format(urllib.parse.quote(structure), type)


def search_np_atlas_by_inchikey_smiles(row):
    # try inchikey - otherwise smiles
    result = get_json_response(np_atlas_url_by_inchikey(row["inchikey"]), True)
    if result is None:
        result = get_json_response(np_atlas_url_by_smiles(row["smiles"]), True)
    if result is None:
        result = get_json_response(np_atlas_url_by_smiles(row["canonical_smiles"]), True)
    # similar structure
    if result is None:
        result = get_json_response(np_atlas_url_similar_structure(row["inchikey"]), True)
    if result is None:
        result = get_json_response(np_atlas_url_similar_structure(row["canonical_smiles"], type="smiles"), True)

    return result


def search_np_atlas(df) -> pd.DataFrame:
    logging.info("Running npatlas search")
    npa = df.apply(lambda row: search_np_atlas_by_inchikey_smiles(row), axis=1)
    df["num_np_atlas_entries"] = npa.apply(len)
    npa = npa.apply(lambda result: result[0] if result else None)

    suffix = NP_ATLAS_SUFFIX
    json_col(df, npa, suffix, "npaid")
    json_col(df, npa, suffix, "original_name")
    json_col(df, npa, suffix, "cluster_id")
    json_col(df, npa, suffix, "node_id")
    json_col(df, npa, suffix, "original_type")
    json_col(df, npa, suffix, "original_organism")
    json_col(df, npa, suffix, "original_doi")

    return df
