import requests
import urllib.parse
import logging

import pandas as pd

from rest_utils import get_json_response_with_headers, json_col

UNI_CHEM_URL = "https://www.ebi.ac.uk/unichem/api/v1/compounds"
UNI_CHEM_SUFFIX = "_unichem"


def search_uni_chem_xref(structure: str, search_type="inchikey"):
    body = {
        "type": search_type,
        "compound": structure
    }
    headers = {"Content-Type": "application/json"}
    result = get_json_response_with_headers(UNI_CHEM_URL, headers, body, True)

    return result


def search_uni_chem_xref_for_row(row):
    # try inchikey - otherwise smiles
    result = search_uni_chem_xref(row["inchikey"], "inchikey")
    return result


def search_uni_chem_xrefs(df) -> pd.DataFrame:
    logging.info("Running uni chem search for cross references")
    npa = df.apply(lambda row: search_uni_chem_xref_for_row(row), axis=1)

    df["uni_chem_compounds"] = [len(result["compounds"]) if result is not None else None
                                for result in npa]

    # npa = [result[0]["sources"] if result is not None and len(result["compounds"]) else None for result in npa]
    df["uni_chem_sources"] = [
        result["compounds"][0]["sources"] if result is not None and len(result["compounds"]) > 0 else None for
        result in npa]

    # suffix = UNI_CHEM_SUFFIX
    # json_col(df, npa, suffix, "npaid")
    # json_col(df, npa, suffix, "original_name")
    # json_col(df, npa, suffix, "cluster_id")
    # json_col(df, npa, suffix, "node_id")
    # json_col(df, npa, suffix, "original_type")
    # json_col(df, npa, suffix, "original_organism")
    # json_col(df, npa, suffix, "original_doi")

    return df
