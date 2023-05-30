import datetime as dt
import logging
# from diskcache import Cache
# cache = Cache("tmpcache")
import pandas as pd
import pubchempy
from joblib import Memory
from pubchempy import Compound, get_compounds

from date_utils import create_expired_entries_dataframe, iso_datetime_now
from meta_constants import MetaColumns
from pandas_utils import notnull, isnull, update_dataframes, get_or_else, make_str_floor_to_int_number
from tqdm import tqdm
from rdkit_mol_identifiers import ensure_smiles_column

tqdm.pandas()

logging.getLogger('pubchempy').setLevel(logging.DEBUG)

memory = Memory("memcache")


def copy_pubchem_synonyms(df: pd.DataFrame) -> pd.DataFrame:
    try:
        pc_synonyms = [pubchem_get_synonyms(compound) for compound in df["pubchem"]]
        df[MetaColumns.synonyms] = [[] if isnull(synonyms) else synonyms for synonyms in df["synonyms"]]
        df[MetaColumns.synonyms] = [a + b for a, b in zip(pc_synonyms, df["synonyms"])]
    except:
        logging.exception("No synonyms")
    return df


def pubchem_search_by_cid(row):
    if notnull(row["pubchem"]):
        return row["pubchem"]

    compound = None
    if "pubchem_cid_parent" in row:
        compound = pubchem_by_cid(row["pubchem_cid_parent"])
    if isnull(compound) and MetaColumns.pubchem_cid in row:
        compound = pubchem_by_cid(row[MetaColumns.pubchem_cid])

    return compound


name_search_columns = [MetaColumns.compound_name, MetaColumns.cas, MetaColumns.input_name, "product_name"]


def pubchem_search_by_names(row) -> str | None:
    if notnull(row["pubchem"]):
        return row["pubchem"]

    for col in name_search_columns:
        value = get_or_else(row, col)
        if isinstance(value, str) and len(value) > 0:
            compound = search_pubchem_by_name(value)
            if notnull(compound):
                return compound

    return None


def transform_pubchem_columns(filtered: pd.DataFrame, apply_structures: bool) -> pd.DataFrame:
    from synonyms import get_first_synonym
    filtered[MetaColumns.pubchem_cid] = [compound.cid for compound in filtered["pubchem"]]
    filtered = make_str_floor_to_int_number(filtered, MetaColumns.pubchem_cid)
    filtered[MetaColumns.pubchem_cid_parent] = filtered[MetaColumns.pubchem_cid]
    filtered[MetaColumns.compound_name] = [get_first_synonym(compound) for compound in filtered["pubchem"]]
    filtered[MetaColumns.iupac] = [compound.iupac_name for compound in filtered["pubchem"]]
    filtered["pubchem_logp"] = [compound.xlogp for compound in filtered["pubchem"]]
    filtered = copy_pubchem_synonyms(filtered)
    if apply_structures:
        filtered[MetaColumns.isomeric_smiles] = [compound.isomeric_smiles for compound in filtered["pubchem"]]
        filtered[MetaColumns.canonical_smiles] = [compound.canonical_smiles for compound in filtered["pubchem"]]

    return filtered


def split_label_structure_sources(df: pd.DataFrame, source_name: str) -> pd.DataFrame:
    input_smiles_df = df[df[MetaColumns.smiles].notnull()].copy()

    if len(input_smiles_df) > 0:
        # concat the new structures with the old ones
        pubchem_smiles_df = df[df[MetaColumns.isomeric_smiles].notnull() |
                               df[MetaColumns.canonical_smiles].notnull()].copy()
        pubchem_smiles_df[MetaColumns.structure_source] = source_name
        pubchem_smiles_df[MetaColumns.smiles] = None
        pubchem_smiles_df = ensure_smiles_column(pubchem_smiles_df)
        # priority pubchem over input later when drop duplicates
        df = pd.concat([pubchem_smiles_df, input_smiles_df], ignore_index=False, sort=False)
    else:
        # no structures were available before - just copy the structures over
        df = ensure_smiles_column(df)
        df[MetaColumns.structure_source] = source_name
    return df


def pubchem_search_structure_by_cid(df: pd.DataFrame, apply_structures: bool,
                                    refresh_expired_entries_after: dt.timedelta = dt.timedelta(
                                        days=90)) -> pd.DataFrame:
    logging.info("Search PubChem by pubchem_cid")
    df["pubchem"] = None

    # only work on expired elements
    # define which rows are old or were not searched before
    filtered = create_expired_entries_dataframe(df, MetaColumns.date_pubchem_cid_search, refresh_expired_entries_after)

    if len(filtered) == 0:  # no need to update
        return df

    # some are filled from the name  or cid search
    filtered["pubchem"] = filtered.progress_apply(lambda row: pubchem_search_by_cid(row), axis=1)
    filtered = filtered[filtered["pubchem"].notnull()].copy()
    # refresh date
    filtered[MetaColumns.date_pubchem_cid_search] = iso_datetime_now()

    # transform create columns, do not copy structures as they are already cleaned by this script
    filtered = transform_pubchem_columns(filtered, apply_structures=apply_structures)

    if apply_structures:
        filtered = split_label_structure_sources(filtered, source_name="pubchem_cid")

    # combine new data with old rows that were not processed
    # keep pubchem to limit name search
    return update_dataframes(filtered, df)
    # return update_dataframes(filtered, df).drop(columns=["pubchem"], errors="ignore")


def pubchem_search_structure_by_name(df: pd.DataFrame, refresh_expired_entries_after: dt.timedelta = dt.timedelta(
    days=90)) -> pd.DataFrame:
    logging.info("Search PubChem by name")
    if "pubchem" not in df.columns:
        df["pubchem"] = None

    # only work on expired elements
    # define which rows are old or were not searched before
    filtered = create_expired_entries_dataframe(df, MetaColumns.date_pubchem_name_search, refresh_expired_entries_after)

    if len(filtered) == 0:  # no need to update
        return df

    # apply search but limit to those without pubchem results from CID search
    filtered = filtered[filtered["pubchem"].isnull()].copy()
    filtered["pubchem"] = filtered.progress_apply(lambda row: pubchem_search_by_names(row), axis=1)
    filtered = filtered[filtered["pubchem"].notnull()].copy()
    # refresh date
    filtered[MetaColumns.date_pubchem_name_search] = iso_datetime_now()

    # transform create new columns
    filtered = transform_pubchem_columns(filtered, apply_structures=True)

    filtered = split_label_structure_sources(filtered, source_name="pubchem_name")

    # combine new data with old rows that were not processed
    # clear pubchem to allow CID and structure search on all of them
    # return update_dataframes(filtered, df)
    return update_dataframes(filtered, df).drop(columns=["pubchem"], errors="ignore")


def pubchem_search_by_structure(df: pd.DataFrame,
                                refresh_expired_entries_after: dt.timedelta = dt.timedelta(days=90)) -> pd.DataFrame:
    logging.info("Search PubChem by structure")
    if "pubchem" not in df.columns:
        df["pubchem"] = None

    # only work on expired elements
    # define which rows are old or were not searched before
    filtered = create_expired_entries_dataframe(df, MetaColumns.date_pubchem_structure_search,
                                                refresh_expired_entries_after)
    if len(filtered) == 0:  # no need to update
        return df

    filtered = filtered[filtered["pubchem"].isnull()].copy()
    filtered["pubchem"] = [search_pubchem_by_structure(smiles, inchi, inchikey) for
                           inchikey, smiles, inchi in
                           zip(filtered["inchikey"], filtered[MetaColumns.smiles], filtered["inchi"])]

    filtered = filtered[filtered["pubchem"].notnull()].copy()
    # refresh date
    filtered[MetaColumns.date_pubchem_structure_search] = iso_datetime_now()

    # transform create columns, do not copy structures as they are already cleaned by this script
    filtered = transform_pubchem_columns(filtered, apply_structures=False)

    # combine new data with old rows that were not processed
    # clear pubchem to allow CID and structure search on all of them
    return update_dataframes(filtered, df).drop(columns=["pubchem"], errors="ignore")


def pubchem_compound_score(comp: Compound):
    """
    SMILES might contain . to denote salts - rate higher those smiles with less . characters
    :param comp:
    :return:
    """
    smiles = comp.canonical_smiles
    if not smiles:
        return 0
    return 1000 - str(smiles).count(".")


@memory.cache
def search_pubchem_by_name(name_or_cas: str) -> Compound | None:
    """
    In pubchem many entries contain the cas as an alternative name - so searching for cas in name works often

    :param name_or_cas: input name or cas
    :return: first compound or None
    """
    if name_or_cas == "NaN":
        return None
    compounds = get_compounds(name_or_cas, "name")
    if not compounds:
        logging.info("Pubchem has no entry named:{}".format(name_or_cas))
        return None
    else:
        compounds.sort(key=lambda comp: pubchem_compound_score(comp), reverse=True)
        return compounds[0]


@memory.cache
def pubchem_by_cid(cid) -> Compound | None:
    return _pubchem_by_cid(cid)


def _pubchem_by_cid(cid, ntry=1, max_tries=10) -> Compound | None:
    try:
        return pubchempy.Compound.from_cid(cid)
    except:
        if ntry < max_tries:
            return _pubchem_by_cid(cid, ntry=ntry + 1, max_tries=max_tries)
        else:
            logging.exception("Cannot retrieve pubchem CID {}".format(cid))
            return None


def search_pubchem_by_structure(smiles=None, inchi=None, inchikey=None) -> Compound | None:
    """
    In pubchem many entries contain the cas as an alternative name - so searching for cas in name works often

    :param smiles:
    :param inchi:
    :param inchikey:
    :return: first compound or None
    """
    if not (smiles or inchi or inchikey):
        return None

    compounds = None
    try:
        if inchikey:
            compounds = get_pubchem_compound(inchikey, "inchikey")
        if not compounds and smiles:
            compounds = get_pubchem_compound(smiles, MetaColumns.smiles)
        if not compounds and inchi:
            compounds = get_pubchem_compound(inchi, "inchi")
        if not compounds:
            logging.info("NO PUBCHEM FOR: smiles:{}  inchi:{}   inchikey:{}".format(smiles, inchi, inchikey))
            return None
        else:
            compounds.sort(key=lambda comp: pubchem_compound_score(comp), reverse=True)
            return compounds[0]
    except:
        logging.warning("FAILED PUBCHEM FOR: smiles:{}  inchi:{}   inchikey:{}".format(smiles, inchi, inchikey))
        return None


@memory.cache
def get_pubchem_compound(value, key):
    return _get_pubchem_compound(value, key)


def _get_pubchem_compound(value, key, ntry=1, max_tries=10):
    try:
        return get_compounds(value, key)
    except:
        if ntry < max_tries:
            return _get_pubchem_compound(value, key, ntry + 1, max_tries)
        else:
            logging.warning("FAILED PUBCHEM FOR: {} (as {})".format(value, key))
            return None


@memory.cache
def pubchem_get_synonyms(compound):
    return _pubchem_get_synonyms(compound)


def _pubchem_get_synonyms(compound, try_n=1, max_tries=10):
    """
    Try to get synonyms with a maximum retry
    :param compound:
    :param try_n: current call
    :param max_tries: maximum tries
    :return: the synonyms or an empty list on maximum number of tries with fail
    """
    if isnull(compound):
        return []
    try:
        synonyms = compound.synonyms
        if synonyms is not None:
            return synonyms
        else:
            return []
    except:
        if try_n < max_tries:
            return _pubchem_get_synonyms(compound, try_n=try_n + 1, max_tries=max_tries)
        else:
            logging.exception("Failed to retrieve synonyms for compound {}".format(compound.cid))
            return []
