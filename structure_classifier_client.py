import logging
import urllib.parse

import pandas as pd
import datetime as dt
from date_utils import create_expired_entries_dataframe, iso_datetime_now
from meta_constants import MetaColumns
from pandas_utils import get_unique_dict, isnull, update_dataframes

# NP_CLASSIFIER_URL = "https://npclassifier.ucsd.edu/classify?smiles={}"
# CLASSYFIRE_URL = "https://gnps-structure.ucsd.edu/classyfire?smiles={}"
from rest_utils import get_json_response, json_col, extract_name, extract_names_array, extract_external_descriptors, \
    join

NP_CLASSIFIER_URL = "https://npclassifier.gnps2.org/classify?smiles={}"
CLASSYFIRE_SMILES_URL = "https://gnps-structure.ucsd.edu/classyfire?smiles={}"
CLASSYFIRE_INCHI_URL = "https://gnps-structure.ucsd.edu/classyfire?inchi={}"
CLASSYFIRE_INCHIKEY_URL = "https://gnps-classyfire.ucsd.edu/entities/{}.json"

CLASSYFIRE_PREFIX = "classyfire"
NP_CLASSIFIER_PREFIX = "npclassifier"


def classyfire_smiles_url(smiles):
    if isnull(smiles):
        return None
    return CLASSYFIRE_SMILES_URL.format(urllib.parse.quote(smiles))


def classyfire_inchikey_url(inchikey):
    if isnull(inchikey):
        return None
    return CLASSYFIRE_INCHIKEY_URL.format(urllib.parse.quote(inchikey))


def np_classifier_url(smiles):
    if isnull(smiles):
        return None
    return NP_CLASSIFIER_URL.format(urllib.parse.quote(smiles))


def apply_np_classifier(df: pd.DataFrame,
                        refresh_expired_entries_after: dt.timedelta = dt.timedelta(days=90)) -> pd.DataFrame:
    """
    :return: a dataframe with the same index as input df, only new columns
    """
    logging.info("Applying np classifier on dataframe")

    # only work on expired elements
    # define which rows are old or were not searched before
    filtered = create_expired_entries_dataframe(df, MetaColumns.date_npclassifier, refresh_expired_entries_after)

    unique_smiles_dict = get_unique_dict(filtered, "smiles")
    for smiles in unique_smiles_dict:
        unique_smiles_dict[smiles] = get_json_response(np_classifier_url(smiles))

    # temp column with json results

    filtered["result_column"] = [unique_smiles_dict[smiles] for smiles in filtered["smiles"]]
    # extract and join values from json array - only isglycoside is already a value
    json_col(filtered, filtered["result_column"], NP_CLASSIFIER_PREFIX, "class_results", join)
    json_col(filtered, filtered["result_column"], NP_CLASSIFIER_PREFIX, "superclass_results", join)
    json_col(filtered, filtered["result_column"], NP_CLASSIFIER_PREFIX, "pathway_results", join)
    json_col(filtered, filtered["result_column"], NP_CLASSIFIER_PREFIX, "isglycoside")
    # json_col(filtered, result_column, NP_CLASSIFIER_PREFIX, "fp1")
    # json_col(filtered, result_column, NP_CLASSIFIER_PREFIX, "fp2")

    # refresh date
    filtered.loc[filtered["result_column"].notnull(), MetaColumns.date_npclassifier] = iso_datetime_now()

    # only keep specific columns
    filter_col = [col for col in filtered.columns if NP_CLASSIFIER_PREFIX in col]
    filtered = filtered[filter_col]
    return filtered
    # return update_dataframes(filtered, df)


def query_classyfire(smiles, inchikey, smiles_dict: dict, inchikey_dict: dict) -> str | None:
    classy_json = smiles_dict.get(smiles)
    if classy_json is not None:
        return classy_json
    classy_json = inchikey_dict.get(inchikey)
    if classy_json is not None:
        return classy_json

    classy_json = get_json_response(classyfire_inchikey_url(inchikey))
    if classy_json is None:
        classy_json = get_json_response(classyfire_smiles_url(smiles))

    smiles_dict[smiles] = classy_json
    inchikey_dict[inchikey] = classy_json
    return classy_json


def apply_classyfire(df: pd.DataFrame,
                     refresh_expired_entries_after: dt.timedelta = dt.timedelta(days=90)) -> pd.DataFrame:
    """
    :return: a dataframe with the same index as input df, only new columns
    """
    logging.info("Applying classyfire on dataframe")

    # only work on expired elements
    # define which rows are old or were not searched before
    filtered = create_expired_entries_dataframe(df, MetaColumns.date_classyfire, refresh_expired_entries_after)

    inchikey_dict = {}
    smiles_dict = {}

    # temp column with json results
    filtered["result_column"] = [query_classyfire(smiles, inchikey, smiles_dict, inchikey_dict) for smiles, inchikey in
                                 zip(filtered["smiles"], filtered["inchikey"])]
    result_column = filtered["result_column"]

    # extract information
    json_col(filtered, result_column, CLASSYFIRE_PREFIX, "kingdom", extract_name)
    json_col(filtered, result_column, CLASSYFIRE_PREFIX, "superclass", extract_name)
    json_col(filtered, result_column, CLASSYFIRE_PREFIX, "class", extract_name)
    json_col(filtered, result_column, CLASSYFIRE_PREFIX, "subclass", extract_name)
    json_col(filtered, result_column, CLASSYFIRE_PREFIX, "intermediate_nodes", extract_names_array)
    json_col(filtered, result_column, CLASSYFIRE_PREFIX, "alternative_parents", extract_names_array)
    json_col(filtered, result_column, CLASSYFIRE_PREFIX, "direct_parent", extract_name)
    json_col(filtered, result_column, CLASSYFIRE_PREFIX, "molecular_framework")
    json_col(filtered, result_column, CLASSYFIRE_PREFIX, "substituents", join)
    json_col(filtered, result_column, CLASSYFIRE_PREFIX, "description")
    json_col(filtered, result_column, CLASSYFIRE_PREFIX, "external_descriptors", extract_external_descriptors)
    json_col(filtered, result_column, CLASSYFIRE_PREFIX, "ancestors", join)
    json_col(filtered, result_column, CLASSYFIRE_PREFIX, "predicted_chebi_terms", join)
    json_col(filtered, result_column, CLASSYFIRE_PREFIX, "predicted_lipidmaps_terms", join)
    json_col(filtered, result_column, CLASSYFIRE_PREFIX, "classification_version")

    # refresh date
    filtered.loc[result_column.notnull(), MetaColumns.date_classyfire] = iso_datetime_now()

    # only keep specific columns
    filter_col = [col for col in filtered.columns if CLASSYFIRE_PREFIX in col]
    filtered = filtered[filter_col]
    return filtered
    # return update_dataframes(filtered, df)
