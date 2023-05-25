import logging
import urllib.parse

import pandas as pd

from pandas_utils import get_unique_dict, isnull

# NP_CLASSIFIER_URL = "https://npclassifier.ucsd.edu/classify?smiles={}"
# CLASSYFIRE_URL = "https://gnps-structure.ucsd.edu/classyfire?smiles={}"
from rest_utils import get_json_response, json_col, extract_name, extract_names_array, extract_external_descriptors, \
    join

NP_CLASSIFIER_URL = "https://npclassifier.gnps2.org/classify?smiles={}"
CLASSYFIRE_SMILES_URL = "https://gnps-structure.ucsd.edu/classyfire?smiles={}"
CLASSYFIRE_INCHI_URL = "https://gnps-structure.ucsd.edu/classyfire?inchi={}"
CLASSYFIRE_INCHIKEY_URL = "https://gnps-classyfire.ucsd.edu/entities/{}.json"

CLASSYFIRE_SUFFIX = "_classyfire"
NP_CLASSIFIER_SUFFIX = "_npclassifier"


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


def apply_np_classifier(df) -> pd.DataFrame:
    logging.info("Applying np classifier on dataframe")
    unique_smiles_dict = get_unique_dict(df, "smiles")
    for smiles in unique_smiles_dict:
        unique_smiles_dict[smiles] = get_json_response(np_classifier_url(smiles))

    # temp column with json results
    result_column = [unique_smiles_dict[smiles] for smiles in df["smiles"]]
    # extract and join values from json array - only isglycoside is already a value
    json_col(df, result_column, NP_CLASSIFIER_SUFFIX, "class_results", join)
    json_col(df, result_column, NP_CLASSIFIER_SUFFIX, "superclass_results", join)
    json_col(df, result_column, NP_CLASSIFIER_SUFFIX, "pathway_results", join)
    json_col(df, result_column, NP_CLASSIFIER_SUFFIX, "isglycoside")
    # json_col(df, result_column, NP_CLASSIFIER_SUFFIX, "fp1")
    # json_col(df, result_column, NP_CLASSIFIER_SUFFIX, "fp2")
    return df


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


def apply_classyfire(df) -> pd.DataFrame:
    logging.info("Applying classyfire on dataframe")
    inchikey_dict = {}
    smiles_dict = {}

    # temp column with json results
    result_column = [query_classyfire(smiles, inchikey, smiles_dict, inchikey_dict) for smiles, inchikey in
                     zip(df["smiles"], df["inchikey"])]

    # extract information
    json_col(df, result_column, CLASSYFIRE_SUFFIX, "kingdom", extract_name)
    json_col(df, result_column, CLASSYFIRE_SUFFIX, "superclass", extract_name)
    json_col(df, result_column, CLASSYFIRE_SUFFIX, "class", extract_name)
    json_col(df, result_column, CLASSYFIRE_SUFFIX, "subclass", extract_name)
    json_col(df, result_column, CLASSYFIRE_SUFFIX, "intermediate_nodes", extract_names_array)
    json_col(df, result_column, CLASSYFIRE_SUFFIX, "alternative_parents", extract_names_array)
    json_col(df, result_column, CLASSYFIRE_SUFFIX, "direct_parent", extract_name)
    json_col(df, result_column, CLASSYFIRE_SUFFIX, "molecular_framework")
    json_col(df, result_column, CLASSYFIRE_SUFFIX, "substituents", join)
    json_col(df, result_column, CLASSYFIRE_SUFFIX, "description")
    json_col(df, result_column, CLASSYFIRE_SUFFIX, "external_descriptors", extract_external_descriptors)
    json_col(df, result_column, CLASSYFIRE_SUFFIX, "ancestors", join)
    json_col(df, result_column, CLASSYFIRE_SUFFIX, "predicted_chebi_terms", join)
    json_col(df, result_column, CLASSYFIRE_SUFFIX, "predicted_lipidmaps_terms", join)
    json_col(df, result_column, CLASSYFIRE_SUFFIX, "classification_version")

    return df
