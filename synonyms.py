import logging
import re
import ast

import pandas as pd
from pyarrow.util import _is_iterable

import pandas_utils as pu
from meta_constants import MetaColumns
from pandas_utils import (
    get_or_else,
    get_unique_list,
    notnull,
    isnull,
    get_first_or_else,
    is_iterable,
)


def ensure_synonyms_column(df: pd.DataFrame) -> pd.DataFrame:
    df["synonyms"] = df.apply(lambda row: get_all_synonyms(row), axis=1)
    return df


def parse_synonyms_to_list(synonyms) -> list:
    if pu.isnull(synonyms):
        return []

    elif isinstance(synonyms, str):
        synonyms = synonyms.strip()
        # parse string from list or from semicolon separated
        if synonyms.startswith("[") and synonyms.endswith("]"):
            return ast.literal_eval(synonyms)

        return [name.strip() for name in synonyms.split(";") if len(name.strip()) > 0]

    elif isinstance(synonyms, list):
        return synonyms
    elif is_iterable(synonyms):
        return [syn for syn in synonyms]


def get_all_synonyms(row):
    synonyms = parse_synonyms_to_list(row["synonyms"] if "synonyms" in row else None)

    new_synonyms = [
        get_or_else(row, "compound_name"),
        get_or_else(row, "input_name"),
        get_or_else(row, "cas"),
    ]

    synonyms += new_synonyms

    synonyms = [x.strip() for x in synonyms if x]
    return get_unique_list(synonyms)


def add_synonyms(old_syn, new_syn, prepend: bool = True) -> list:
    """
    Add synonyms from two lists together
    :param old_syn:
    :param new_syn:
    :param prepend:
    :return: unique list
    """
    if isnull(old_syn):
        old_syn = []
    elif isinstance(old_syn, str):
        old_syn = parse_synonyms_to_list(old_syn)
    if isnull(new_syn):
        new_syn = []
    elif isinstance(new_syn, str):
        new_syn = parse_synonyms_to_list(new_syn)
    comb = new_syn + old_syn if prepend else old_syn + new_syn
    comb = [x.strip() for x in comb if x]
    return get_unique_list(comb)


def add_synonyms_columns(
    df, new_syn_column_header: str = None, new_synonyms=None, prepend: bool = False
) -> pd.DataFrame:
    """
    Add all synonyms to synonyms column of df
    :param df:
    :param new_syn_column_header:
    :return: column will contain unique lists, returns the input df
    """
    if notnull(new_syn_column_header) and new_syn_column_header in df:
        col = df[new_syn_column_header]
        df[MetaColumns.synonyms] = [
            add_synonyms(olds, news, prepend)
            for olds, news in zip(df[MetaColumns.synonyms], col)
        ]

    if notnull(new_synonyms):
        col = new_synonyms
        df[MetaColumns.synonyms] = [
            add_synonyms(olds, news, prepend)
            for olds, news in zip(df[MetaColumns.synonyms], col)
        ]

    return df


def find_unii(synonyms):
    unii_generator = (
        re.sub("[ .;:\-]|UNII", "", name.upper())
        for name in synonyms
        if "UNII" in name.upper()
    )
    return next(unii_generator, None)


def find_schembl(synonyms):
    schembl_generator = (
        name.upper() for name in synonyms if name.upper().startswith("SCHEMBL")
    )
    return next(schembl_generator, None)


def find_chembl_id(synonyms):
    chembl_generator = (
        name.upper() for name in synonyms if name.upper().startswith("CHEMBL")
    )
    return next(chembl_generator, None)


def find_zinc(synonyms):
    zinc_generator = (
        name.upper() for name in synonyms if name.upper().startswith("ZINC")
    )
    return next(zinc_generator, None)


def find_drugbank(synonyms):
    for s in synonyms:
        drug = cleanup_drugbank_id(s)
        if drug:
            return drug
    return None
    # drugbank_generator = (cleanup_drugbank_id(name) for name in synonyms)
    # return next((db_id for db_id in drugbank_generator if db_id), None)


def cleanup_drugbank_id(input):
    pattern = "^DB.*\d"
    anti_pattern = "[ACE-Z]"
    input = input.upper()
    if re.search(pattern, input) and not re.search(anti_pattern, input):
        return re.sub("[^0-9DB]", "", input)
    else:
        return None


def extract_synonym_ids(df: pd.DataFrame) -> pd.DataFrame:
    source = pd.DataFrame()
    source["unii"] = [find_unii(synonyms) for synonyms in df["synonyms"]]
    source["schembl_id"] = [find_schembl(synonyms) for synonyms in df["synonyms"]]
    source["chembl_id"] = [find_chembl_id(synonyms) for synonyms in df["synonyms"]]
    source["zinc_id"] = [find_zinc(synonyms) for synonyms in df["synonyms"]]
    source["drugbank_id"] = [find_drugbank(synonyms) for synonyms in df["synonyms"]]

    return pu.combine_dfs_fill_missing_values(df, source)


def use_first_synonym_as_compound_name(df: pd.DataFrame) -> pd.DataFrame:
    if MetaColumns.synonyms not in df.columns:
        return df

    df[MetaColumns.compound_name] = [
        get_first_or_else(synonyms) for synonyms in df[MetaColumns.synonyms]
    ]
    return df
