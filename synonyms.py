import logging
import re

import pandas as pd

from pandas_utils import get_or_else, get_unique_list
from pubchem_client import pubchem_get_synonyms


def ensure_synonyms_column(df: pd.DataFrame) -> pd.DataFrame:
    df["synonyms"] = df.apply(lambda row: get_all_synonyms(row), axis=1)
    return df


def get_all_synonyms(row):
    synonyms = [
        get_or_else(row, "product_name"),
        get_or_else(row, "cas"),
        get_or_else(row, "compound_name"),
        get_or_else(row, "input_compound_name"),
    ]

    old = row["synonyms"] if "synonyms" in row else None
    try:
        if isinstance(old, str):
            synonyms.append(old)
        elif old is not None:
            synonyms = synonyms + old
    except:
        logging.exception("Cannot concat synonyms")

    # TODO remove this Synonyms columns?
    # synonyms.extend([s.strip() for s in str(get_or_else(row, "synonyms", "")).split(";")])

    synonyms = [x.strip() for x in synonyms if x]
    return get_unique_list(synonyms)


def get_first_synonym(compound):
    synonyms = pubchem_get_synonyms(compound)
    if synonyms is None or len(synonyms) <= 0:
        return None
    return synonyms[0]


def find_unii(synonyms):
    unii_generator = (re.sub('[ .;:\-]|UNII', '', name.upper()) for name in synonyms if "UNII" in name.upper())
    return next(unii_generator, None)


def find_schembl(synonyms):
    schembl_generator = (name.upper() for name in synonyms if name.upper().startswith("SCHEMBL"))
    return next(schembl_generator, None)


def find_chembl_id(synonyms):
    chembl_generator = (name.upper() for name in synonyms if name.upper().startswith("CHEMBL"))
    return next(chembl_generator, None)


def find_zinc(synonyms):
    zinc_generator = (name.upper() for name in synonyms if name.upper().startswith("ZINC"))
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
    df["unii"] = [find_unii(synonyms) for synonyms in df["synonyms"]]
    df["schembl_id"] = [find_schembl(synonyms) for synonyms in df["synonyms"]]
    df["chembl_id"] = [find_chembl_id(synonyms) for synonyms in df["synonyms"]]
    df["zinc_id"] = [find_zinc(synonyms) for synonyms in df["synonyms"]]
    df["drugbank_id"] = [find_drugbank(synonyms) for synonyms in df["synonyms"]]

    return df
