import logging
import urllib.parse

import pandas as pd
import datetime as dt
from date_utils import create_expired_entries_dataframe, iso_datetime_now
from meta_constants import MetaColumns
from pandas_utils import (
    get_unique_dict,
    isnull,
    update_dataframes,
    notnull,
    notnull_not_empty,
    isnull_or_empty,
)
from tqdm import tqdm
from joblib import Memory

memory = Memory("memcache")

# NP_CLASSIFIER_URL = "https://npclassifier.ucsd.edu/classify?smiles={}"
# CLASSYFIRE_URL = "https://gnps-structure.ucsd.edu/classyfire?smiles={}"
from rest_utils import (
    get_json_response,
    json_col,
    extract_name,
    extract_names_array,
    extract_external_descriptors,
    join,
)

NP_CLASSIFIER_URL = "https://npclassifier.gnps2.org/classify?smiles={}"
CLASSYFIRE_SMILES_URL = "https://structure.gnps2.org/classyfire?smiles={}"
CLASSYFIRE_INCHI_URL = "https://structure.gnps2.org/classyfire?inchi={}"
CLASSYFIRE_INCHIKEY_URL = "https://classyfire.gnps2.org/entities/{}.json"

CLASSYFIRE_PREFIX = "classyfire"
NP_CLASSIFIER_PREFIX = "npclassifier"


def classyfire_smiles_url(smiles):
    if isnull_or_empty(smiles):
        return None
    return CLASSYFIRE_SMILES_URL.format(urllib.parse.quote(smiles))


def classyfire_inchikey_url(inchikey):
    if isnull_or_empty(inchikey):
        return None
    return CLASSYFIRE_INCHIKEY_URL.format(urllib.parse.quote(inchikey))


def np_classifier_url(smiles):
    if isnull_or_empty(smiles):
        return None
    return NP_CLASSIFIER_URL.format(urllib.parse.quote(smiles))


def apply_np_classifier(
    df: pd.DataFrame,
    refresh_expired_entries_after: dt.timedelta = dt.timedelta(days=90),
) -> pd.DataFrame:
    """
    :return: a dataframe with the same index as input df, only new columns
    """
    logging.info("Applying np classifier on dataframe")

    # only work on expired elements
    # define which rows are old or were not searched before
    filtered = create_expired_entries_dataframe(
        df, MetaColumns.date_npclassifier, refresh_expired_entries_after
    )

    unique_smiles_dict = get_unique_dict(filtered, MetaColumns.smiles)
    for smiles in tqdm(unique_smiles_dict):
        unique_smiles_dict[smiles] = query_np_classifier(smiles)

    # temp column with json results

    filtered["result_column"] = [
        unique_smiles_dict[smiles] for smiles in filtered[MetaColumns.smiles]
    ]
    # extract and join values from json array - only isglycoside is already a value
    json_col(
        filtered, filtered["result_column"], NP_CLASSIFIER_PREFIX, "class_results", join
    )
    json_col(
        filtered,
        filtered["result_column"],
        NP_CLASSIFIER_PREFIX,
        "superclass_results",
        join,
    )
    json_col(
        filtered,
        filtered["result_column"],
        NP_CLASSIFIER_PREFIX,
        "pathway_results",
        join,
    )
    json_col(filtered, filtered["result_column"], NP_CLASSIFIER_PREFIX, "isglycoside")
    # json_col(filtered, result_column, NP_CLASSIFIER_PREFIX, "fp1")
    # json_col(filtered, result_column, NP_CLASSIFIER_PREFIX, "fp2")

    # refresh date
    filtered.loc[
        filtered["result_column"].notnull(), MetaColumns.date_npclassifier
    ] = iso_datetime_now()

    # only keep specific columns
    filter_col = [col for col in filtered.columns if NP_CLASSIFIER_PREFIX in col]
    filtered = filtered[filter_col]
    return filtered
    # return update_dataframes(filtered, df)


def query_np_classifier(smiles) -> str | None:
    try:
        return _query_np_classifier(smiles)
    except:
        logging.info(f"Failed np classifier for {smiles}")
        return None


@memory.cache
def _query_np_classifier(smiles) -> str | None:
    if isnull_or_empty(smiles):
        return None

    response = get_json_response(np_classifier_url(smiles))
    if response is None:
        raise Exception()
    return response


def query_classyfire(smiles, inchikey) -> str | None:
    try:
        return _query_classyfire(smiles, inchikey)
    except:
        logging.info(f"Failed classyfire for {smiles} and {inchikey}")
        return None


@memory.cache
def _query_classyfire(smiles, inchikey) -> str | None:
    if isnull_or_empty(smiles) and isnull_or_empty(inchikey):
        return None

    # TODO reduce timeout again to 2 and 3 sec
    classy_json = None
    if notnull_not_empty(inchikey):
        classy_json = get_json_response(classyfire_inchikey_url(inchikey), timeout=5)

    if classy_json is None and notnull_not_empty(smiles):
        classy_json = get_json_response(classyfire_smiles_url(smiles), timeout=5)

    if classy_json is None:
        raise Exception()
    return classy_json


def apply_classyfire(
    df: pd.DataFrame,
    refresh_expired_entries_after: dt.timedelta = dt.timedelta(days=90),
) -> pd.DataFrame:
    """
    :return: a dataframe with the same index as input df, only new columns
    """
    logging.info("Applying classyfire on dataframe")

    # only work on expired elements
    # define which rows are old or were not searched before
    filtered = create_expired_entries_dataframe(
        df, MetaColumns.date_classyfire, refresh_expired_entries_after
    )

    # temp column with json results
    filtered["result_column"] = [
        query_classyfire(smiles, inchikey)
        for smiles, inchikey in tqdm(
            zip(filtered["smiles"], filtered["inchikey"]), total=len(filtered)
        )
    ]
    result_column = filtered["result_column"]

    # extract information
    json_col(filtered, result_column, CLASSYFIRE_PREFIX, "kingdom", extract_name)
    json_col(filtered, result_column, CLASSYFIRE_PREFIX, "superclass", extract_name)
    json_col(filtered, result_column, CLASSYFIRE_PREFIX, "class", extract_name)
    json_col(filtered, result_column, CLASSYFIRE_PREFIX, "subclass", extract_name)
    json_col(
        filtered,
        result_column,
        CLASSYFIRE_PREFIX,
        "intermediate_nodes",
        extract_names_array,
    )
    json_col(
        filtered,
        result_column,
        CLASSYFIRE_PREFIX,
        "alternative_parents",
        extract_names_array,
    )
    json_col(filtered, result_column, CLASSYFIRE_PREFIX, "direct_parent", extract_name)
    json_col(filtered, result_column, CLASSYFIRE_PREFIX, "molecular_framework")
    json_col(filtered, result_column, CLASSYFIRE_PREFIX, "substituents", join)
    json_col(filtered, result_column, CLASSYFIRE_PREFIX, "description")
    json_col(
        filtered,
        result_column,
        CLASSYFIRE_PREFIX,
        "external_descriptors",
        extract_external_descriptors,
    )
    json_col(filtered, result_column, CLASSYFIRE_PREFIX, "ancestors", join)
    json_col(filtered, result_column, CLASSYFIRE_PREFIX, "predicted_chebi_terms", join)
    json_col(
        filtered, result_column, CLASSYFIRE_PREFIX, "predicted_lipidmaps_terms", join
    )
    json_col(filtered, result_column, CLASSYFIRE_PREFIX, "classification_version")

    # refresh date
    filtered.loc[
        result_column.notnull(), MetaColumns.date_classyfire
    ] = iso_datetime_now()

    # only keep specific columns
    filter_col = [col for col in filtered.columns if CLASSYFIRE_PREFIX in col]
    filtered = filtered[filter_col]
    return filtered
    # return update_dataframes(filtered, df)
