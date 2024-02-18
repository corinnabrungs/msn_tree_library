import logging
import urllib.parse
import pandas as pd

import pandas_utils
from date_utils import iso_datetime_now, create_expired_entries_dataframe
from meta_constants import MetaColumns
from pandas_utils import remove_empty_lists_values, isnull, isnull_or_empty
from rest_utils import get_json_response, json_col
import datetime as dt
from joblib import Memory

memory = Memory("memcache")

# NPAtlas:
# ORIGIN ORGANISM TYPE
NP_ATLAS_URL = "https://www.npatlas.org/api/v1/compounds/basicSearch?method=full&{}&threshold=0&orderby=npaid&ascending=true&limit=10"
NP_ATLAS_STRUCTURE_SEARCH_URL = (
    "https://www.npatlas.org/api/v1/compounds/structureSearch?structure={}"
    "&type={}&method=sim&threshold=0.9999&skip=0&limit=10&stereo=false"
)
NP_ATLAS_PREFIX = "npatlas"


def np_atlas_url_by_exact_structure(structure, search_type="inchikey"):
    return NP_ATLAS_URL.format(search_type + "=" + urllib.parse.quote(structure))


def np_atlas_url_similar_structure(structure, search_type="inchikey"):
    return NP_ATLAS_URL.format(urllib.parse.quote(structure), search_type)


def npatlas_url(
    structure: str, search_type: str = "inchikey", similarity: str = "exact"
):
    url = None
    if similarity == "exact":
        url = np_atlas_url_by_exact_structure(structure, search_type)
    if similarity == "similar":
        url = np_atlas_url_similar_structure(structure, search_type)
    if url is None:
        raise ValueError(
            "url is None this means npatlas search had wrong search or similarity type"
        )
    return url


def search_np_atlas_by_inchikey_smiles(row, search_similar: bool = False):
    # try inchikey - otherwise smiles
    similarity = "similar" if search_similar else "exact"
    result = None
    if "inchikey" in row:
        result = query_npatlas(
            row["inchikey"], search_type="inchikey", similarity=similarity
        )
    if result is None and "smiles" in row:
        result = query_npatlas(
            row["smiles"], search_type="smiles", similarity=similarity
        )
    if result is None and "canonical_smiles" in row:
        result = query_npatlas(
            row["canonical_smiles"], search_type="smiles", similarity=similarity
        )

    return result


def query_npatlas(
    structure: str, search_type: str = "inchikey", similarity: str = "exact"
):
    try:
        return _query_npatlas(structure, search_type, similarity)
    except:
        logging.info(
            f"Failed npatlas for {structure} as {search_type} with sim {similarity}"
        )
        return None


@memory.cache
def _query_npatlas(
    structure: str, search_type: str = "inchikey", similarity: str = "exact"
):
    """
    :param structure: input structure
    :param search_type: inchikey or smiles
    :param similarity: exact or similar
    :return: list of matches or None if failed
    """
    if isnull_or_empty(structure):
        return None

    url = npatlas_url(structure, search_type, similarity)

    result = get_json_response(url, True, timeout=3)
    if result is None:
        raise Exception("Failed npatlas service")
    return result


def search_np_atlas(
    df: pd.DataFrame,
    refresh_expired_entries_after: dt.timedelta = dt.timedelta(days=90),
) -> pd.DataFrame:
    logging.info("Running npatlas search")

    # only work on expired elements
    # define which rows are old or were not searched before
    filtered = create_expired_entries_dataframe(
        df, MetaColumns.date_npatlas, refresh_expired_entries_after
    )

    filtered["results"] = filtered.apply(
        lambda row: search_np_atlas_by_inchikey_smiles(row), axis=1
    )
    # refresh date
    filtered.loc[
        filtered["results"].notnull(), MetaColumns.date_npatlas
    ] = iso_datetime_now()

    filtered = remove_empty_lists_values(filtered, "results")  # empty results
    filtered = filtered[filtered["results"].notnull()].copy()
    npa = filtered["results"]

    prefix = NP_ATLAS_PREFIX
    filtered[prefix + "_num_entries"] = npa.apply(len)
    npa = npa.apply(lambda result: result[0])

    json_col(filtered, npa, prefix, "npaid")
    json_col(filtered, npa, prefix, "original_name")
    json_col(filtered, npa, prefix, "cluster_id")
    json_col(filtered, npa, prefix, "node_id")
    json_col(filtered, npa, prefix, "original_type")
    json_col(filtered, npa, prefix, "original_organism")
    json_col(filtered, npa, prefix, "original_doi")

    # only keep specific columns
    filter_col = [col for col in filtered.columns if prefix in col]
    return filtered[filter_col]
