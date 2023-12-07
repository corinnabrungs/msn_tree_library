import datetime as dt
import logging

# from diskcache import Cache
# cache = Cache("tmpcache")
import pandas as pd
import pubchempy
from joblib import Memory
from pubchempy import Compound, get_compounds

import pandas_utils
import synonyms

from date_utils import create_expired_entries_dataframe, iso_datetime_now
from meta_constants import MetaColumns
from pandas_utils import (
    notnull,
    isnull,
    update_dataframes,
    get_or_else,
    make_str_floor_to_int_number,
    get_first_or_else,
)
from tqdm import tqdm
from rdkit_mol_identifiers import ensure_smiles_column

tqdm.pandas()

logging.getLogger("pubchempy").setLevel(logging.DEBUG)

memory = Memory("memcache")


def prepend_pubchem_synonyms(df: pd.DataFrame) -> pd.DataFrame:
    try:
        pc_synonyms = [pubchem_get_synonyms(compound) for compound in df["pubchem"]]
        df = synonyms.ensure_synonyms_column(df)
        df[MetaColumns.synonyms] = [
            synonyms.add_synonyms(old, new)
            for new, old in zip(pc_synonyms, df["synonyms"])
        ]
    except:
        logging.exception("No synonyms")
    return df


cid_search_columns = [MetaColumns.pubchem_cid, MetaColumns.input_pubchem_cid]


def pubchem_search_by_cid(row):
    if notnull(row["pubchem"]):
        return row["pubchem"]

    for col in cid_search_columns:
        if col in row:
            comp = pubchem_by_cid(row[col])
            if notnull(comp):
                return comp

    return None


name_search_columns = [
    MetaColumns.compound_name,
    MetaColumns.cas,
    MetaColumns.input_name,
    "product_name",
]


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


def transform_pubchem_columns(
    filtered: pd.DataFrame, apply_structures: bool
) -> pd.DataFrame:
    filtered[MetaColumns.pubchem_cid] = [
        compound.cid for compound in filtered["pubchem"]
    ]
    filtered = make_str_floor_to_int_number(filtered, MetaColumns.pubchem_cid)
    filtered[MetaColumns.iupac] = [
        compound.iupac_name for compound in filtered["pubchem"]
    ]
    filtered["pubchem_logp"] = [compound.xlogp for compound in filtered["pubchem"]]
    filtered = prepend_pubchem_synonyms(filtered)
    filtered[MetaColumns.compound_name] = [
        get_first_or_else(synonyms) for synonyms in filtered[MetaColumns.synonyms]
    ]
    if apply_structures:
        filtered[MetaColumns.isomeric_smiles] = [
            compound.isomeric_smiles for compound in filtered["pubchem"]
        ]
        filtered[MetaColumns.canonical_smiles] = [
            compound.canonical_smiles for compound in filtered["pubchem"]
        ]

    return filtered


def split_label_structure_sources(df: pd.DataFrame, source_name: str) -> pd.DataFrame:
    input_smiles_df = df[df[MetaColumns.smiles].notnull()].copy()

    if len(input_smiles_df) > 0:
        # concat the new structures with the old ones
        pubchem_smiles_df = df[
            df[MetaColumns.isomeric_smiles].notnull()
            | df[MetaColumns.canonical_smiles].notnull()
        ].copy()
        pubchem_smiles_df[MetaColumns.structure_source] = source_name
        pubchem_smiles_df[MetaColumns.smiles] = None
        pubchem_smiles_df = ensure_smiles_column(pubchem_smiles_df)
        # priority pubchem over input later when drop duplicates
        df = pd.concat(
            [pubchem_smiles_df, input_smiles_df], ignore_index=False, sort=False
        )
    else:
        # no structures were available before - just copy the structures over
        df = ensure_smiles_column(df)
        df[MetaColumns.structure_source] = source_name
    return df


def pubchem_search_structure_by_cid(
    df: pd.DataFrame,
    apply_structures: bool,
    refresh_expired_entries_after: dt.timedelta = dt.timedelta(days=90),
) -> pd.DataFrame:
    logging.info("Search PubChem by pubchem_cid")
    df["pubchem"] = None

    # only work on expired elements
    # define which rows are old or were not searched before
    filtered = create_expired_entries_dataframe(
        df, MetaColumns.date_pubchem_cid_search, refresh_expired_entries_after
    )

    if len(filtered) == 0:  # no need to update
        return df

    # some are filled from the name  or cid search
    filtered["pubchem"] = filtered.progress_apply(
        lambda row: pubchem_search_by_cid(row), axis=1
    )
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


def pubchem_search_parent(
    df: pd.DataFrame,
    apply_structures: bool,
    refresh_expired_entries_after: dt.timedelta = dt.timedelta(days=90),
) -> pd.DataFrame:
    logging.info("Search PubChem parents pubchem_cid")
    df["pubchem"] = None

    # only work on expired elements
    # define which rows are old or were not searched before
    filtered = create_expired_entries_dataframe(
        df, MetaColumns.date_pubchem_parent_cid_search, refresh_expired_entries_after
    )

    if len(filtered) == 0:  # no need to update
        return df

    # some are filled from the name  or cid search
    if MetaColumns.input_pubchem_cid not in filtered:
        filtered[MetaColumns.input_pubchem_cid] = None
    if MetaColumns.pubchem_cid not in filtered:
        filtered[MetaColumns.pubchem_cid] = None
    filtered[MetaColumns.input_pubchem_cid] = [
        input_cid if notnull(input_cid) else cid
        for input_cid, cid in zip(
            filtered[MetaColumns.input_pubchem_cid], filtered[MetaColumns.pubchem_cid]
        )
    ]
    # get parent CIDs, filter out failed, map to compound, apply all columns
    filtered[MetaColumns.pubchem_cid] = filtered[
        MetaColumns.input_pubchem_cid
    ].progress_apply(lambda cid: get_pubchem_parent_cid(cid))
    filtered = filtered[filtered[MetaColumns.pubchem_cid].notnull()].copy()
    filtered["pubchem"] = filtered[MetaColumns.pubchem_cid].progress_apply(
        lambda cid: pubchem_by_cid(cid)
    )
    filtered = filtered[filtered["pubchem"].notnull()].copy()
    # refresh date
    filtered[MetaColumns.date_pubchem_parent_cid_search] = iso_datetime_now()

    # transform create columns, do not copy structures as they are already cleaned by this script
    filtered = transform_pubchem_columns(filtered, apply_structures=apply_structures)

    if apply_structures:
        filtered = split_label_structure_sources(
            filtered, source_name="parent_pubchem_cid"
        )

    # combine new data with old rows that were not processed
    # keep pubchem to limit name search
    return update_dataframes(filtered, df)


def pubchem_search_structure_by_name(
    df: pd.DataFrame,
    refresh_expired_entries_after: dt.timedelta = dt.timedelta(days=90),
) -> pd.DataFrame:
    logging.info("Search PubChem by name")
    if "pubchem" not in df.columns:
        df["pubchem"] = None

    # only work on expired elements
    # define which rows are old or were not searched before
    filtered = create_expired_entries_dataframe(
        df, MetaColumns.date_pubchem_name_search, refresh_expired_entries_after
    )

    if len(filtered) == 0:  # no need to update
        return df

    # apply search but limit to those without pubchem results from CID search
    filtered = filtered[filtered["pubchem"].isnull()].copy()
    filtered["pubchem"] = filtered.progress_apply(
        lambda row: pubchem_search_by_names(row), axis=1
    )
    filtered = filtered[filtered["pubchem"].notnull()].copy()
    # refresh date
    filtered[MetaColumns.date_pubchem_name_search] = iso_datetime_now()

    # transform create new columns
    filtered = transform_pubchem_columns(filtered, apply_structures=True)

    filtered = split_label_structure_sources(filtered, source_name="pubchem_name")

    # combine new data with old rows that were not processed
    # clear pubchem to allow CID and structure search on all of them
    # return update_dataframes(filtered, df)
    df = update_dataframes(filtered, df).drop(columns=["pubchem"], errors="ignore")
    df = make_str_floor_to_int_number(df, MetaColumns.pubchem_cid)
    return df


def pubchem_search_by_structure(
    df: pd.DataFrame,
    refresh_expired_entries_after: dt.timedelta = dt.timedelta(days=90),
) -> pd.DataFrame:
    logging.info("Search PubChem by structure")
    if "pubchem" not in df.columns:
        df["pubchem"] = None

    # only work on expired elements
    # define which rows are old or were not searched before
    filtered = create_expired_entries_dataframe(
        df, MetaColumns.date_pubchem_structure_search, refresh_expired_entries_after
    )
    if len(filtered) == 0:  # no need to update
        return df

    filtered = filtered[filtered["pubchem"].isnull()].copy()
    filtered["pubchem"] = [
        search_pubchem_by_structure(smiles, inchi, inchikey)
        for inchikey, smiles, inchi in zip(
            filtered["inchikey"], filtered[MetaColumns.smiles], filtered["inchi"]
        )
    ]

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
    try:
        compounds = get_compounds(name_or_cas, "name")
        if not compounds:
            logging.info("Pubchem has no entry named:{}".format(name_or_cas))
            return None
        else:
            compounds.sort(key=lambda comp: pubchem_compound_score(comp), reverse=True)
            return compounds[0]
    except:
        logging.warning(f"FAILED PUBCHEM by name {name_or_cas}")
        pass


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
            logging.exception("Cannot retrieve PUBCHEM CID {}".format(cid))
            return None


def search_pubchem_by_structure(
    smiles=None, inchi=None, inchikey=None
) -> Compound | None:
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
            compounds = get_pubchem_compound(inchikey, MetaColumns.inchikey)
        if not compounds and smiles:
            compounds = get_pubchem_compound(smiles, MetaColumns.smiles)
        if not compounds and inchi:
            compounds = get_pubchem_compound(inchi, MetaColumns.inchi)
        if not compounds:
            logging.info(
                "NO PUBCHEM FOR: smiles:{}  inchi:{}   inchikey:{}".format(
                    smiles, inchi, inchikey
                )
            )
            return None
        else:
            compounds.sort(key=lambda comp: pubchem_compound_score(comp), reverse=True)
            return compounds[0]
    except:
        logging.warning(
            "FAILED PUBCHEM FOR: smiles:{}  inchi:{}   inchikey:{}".format(
                smiles, inchi, inchikey
            )
        )
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
            logging.exception(
                "FAILED to retrieve PUBCHEM synonyms for compound {}".format(
                    compound.cid
                )
            )
            return []


@memory.cache
def get_pubchem_parent_cid(cid, orphans_as_self=True) -> str | None:
    return _get_pubchem_parent_cid(cid, orphans_as_self)


def _get_pubchem_parent_cid(
    cid, orphans_as_self=True, try_n=1, max_tries=10
) -> str | None:
    """
    From a pubchem_cid, retreive the parent compound's cid.
    If function is unsuccesful in retrieving a single parent,
    `orphans_as_self = True` returns `cid` rather than None.

    According to pubmed:

    > A parent is conceptually the "important" part of the molecule
    > when the molecule has more than one covalent component.
    > Specifically, a parent component must have at least one carbon
    > and contain at least 70% of the heavy (non-hydrogen) atoms of
    > all the unique covalent units (ignoring stoichiometry).
    > Note that this is a very empirical definition and is subject to change.
    A parallel query can be executed using the REST PUG API:
    http://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/11477084/cids/XML?cids_type=parent
    """
    if isnull(cid):
        return None

    try:
        parent_cids = pubchempy.get_cids(
            identifier=cid, namespace="cid", domain="compound", cids_type="parent"
        )
    except pubchempy.BadRequestError:
        logging.info("Error getting parent of {}".format(cid))
        return None
    except Exception:
        if try_n < max_tries:
            return _get_pubchem_parent_cid(cid, orphans_as_self, try_n + 1, max_tries)
        else:
            logging.exception("Error getting parent of {}".format(cid))
            return None  # error return None to redo later
    try:
        if len(parent_cids) > 0:
            return str(parent_cids[0])
        return cid if orphans_as_self else None
    except Exception:
        logging.exception(
            "Error getting parent of {}. Parents retrieved: {}".format(cid, parent_cids)
        )
        return cid if orphans_as_self else None
