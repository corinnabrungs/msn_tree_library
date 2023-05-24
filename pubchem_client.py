import logging
# from diskcache import Cache
# cache = Cache("tmpcache")
import numpy as np
import pandas as pd
import pubchempy
from joblib import Memory
from pubchempy import Compound, get_compounds

from pandas_utils import notnull, isnull
from tqdm import tqdm

tqdm.pandas()

logging.getLogger('pubchempy').setLevel(logging.DEBUG)

memory = Memory("memcache")


def pubchem_search_by_cid(row):
    if notnull(row["pubchem"]):
        return row["pubchem"]

    compound = None
    if "pubchem_cid_parent" in row:
        compound = pubchem_by_cid(row["pubchem_cid_parent"])
    if isnull(compound) and "pubchem_cid" in row:
        compound = pubchem_by_cid(row["pubchem_cid"])

    return compound


def pubchem_search_by_names(row) -> str | None:
    if notnull(row["pubchem"]):
        return row["pubchem"]

    compound = pubchem_search_by_cid(row)
    if isnull(compound) and "compound_name" in row and notnull(row["compound_name"]):
        compound = search_pubchem_by_name(str(row["compound_name"]))
    if isnull(compound) and "cas" in row and notnull(row["cas"]):
        compound = search_pubchem_by_name(row["cas"])
    if isnull(compound) and "product_name" in row:
        compound = search_pubchem_by_name(row["product_name"])
    # only one compound was found as CAS-
    if isnull(compound) and "cas" in row:
        compound = search_pubchem_by_name("CAS-{}".format(row["product_name"]))

    return compound


def pubchem_search_structure_by_name(df) -> pd.DataFrame:
    from synonyms import get_first_synonym

    logging.info("Search PubChem by name")
    df["pubchem"] = None
    df["pubchem"] = df.progress_apply(lambda row: pubchem_search_by_names(row), axis=1)

    df["pubchem_cid"] = pd.array([compound.cid if pd.notnull(compound) else np.NAN for compound in df["pubchem"]],
                                 dtype=pd.Int64Dtype())
    df["isomerical_smiles"] = [compound.isomeric_smiles if pd.notnull(compound) else np.NAN for compound in
                               df["pubchem"]]
    df["canonical_smiles"] = [compound.canonical_smiles if pd.notnull(compound) else np.NAN for compound in
                              df["pubchem"]]

    df["pubchem_cid_parent"] = df["pubchem_cid"]
    df["compound_name"] = [get_first_synonym(compound) for compound in df["pubchem"]]
    df["iupac"] = [compound.iupac_name if pd.notnull(compound) else np.NAN for compound in df["pubchem"]]
    try:
        df["synonyms"] = df["synonyms"] + [pubchem_get_synonyms(compound) if pd.notnull(compound) else [] for compound
                                           in df["pubchem"]]
    except:
        logging.exception("No synonyms")
    df["pubchem_logp"] = [compound.xlogp if pd.notnull(compound) else np.NAN for compound in df["pubchem"]]

    # TODO remove this part as we are not dropping any columns
    # # drop extra columns
    # columns_to_keep = df.columns.isin(
    #     ["pubchem", "Cat. No.", "product_name", "synonyms", "cas", "smiles", "pubchem_cid", "pubchem_cid_parent",
    #      "isomerical_smiles",
    #      "canonical_smiles", "mixed_location_plate1", "lib_plate_well", "URL", "Target", "Information", "Pathway",
    #      "Research Area", "Clinical Information", "gnps_libid", "compound_name", "iupac", "pubchem_logp", "entries"])
    #
    # df = df[df.columns[columns_to_keep]]
    # concat the new structures with the old ones
    if "smiles" in df and len(df[df["smiles"].notna()]) > 0:
        dfa = df.copy()

        dfa["source"] = "input"
        dfb = df.copy()
        dfb["smiles"] = [iso if notnull(iso) else smiles for iso, smiles in zip(df["isomerical_smiles"], df["smiles"])]
        dfb["source"] = "PubChem"

        dfb = dfb[dfb["smiles"].notna()]

        df = pd.concat([dfb, dfa], ignore_index=True, sort=False)

    else:
        df["smiles"] = [iso if notnull(iso) else smiles for iso, smiles in zip(df["isomerical_smiles"], df["smiles"])]
        df["source"] = "PubChem"
    return df


def pubchem_search_by_structure(df) -> pd.DataFrame:
    from synonyms import get_first_synonym

    logging.info("Search PubChem by structure")
    df["pubchem"] = df.progress_apply(lambda row: pubchem_search_by_cid(row), axis=1)

    df["pubchem"] = [search_pubchem_by_structure(smiles, inchi, inchikey) if isnull(compound) else compound for
                     compound, inchikey, smiles, inchi in
                     zip(df["pubchem"], df["inchikey"], df["smiles"], df["inchi"])]

    df["pubchem_cid_parent"] = pd.array(
        [compound.cid if pd.notnull(compound) else np.NAN for compound in df["pubchem"]],
        dtype=pd.Int64Dtype())
    df["compound_name"] = [get_first_synonym(compound) for compound in df["pubchem"]]
    df["iupac"] = [compound.iupac_name if pd.notnull(compound) else np.NAN for compound in df["pubchem"]]
    try:
        df["synonyms"] = df["synonyms"] + [pubchem_get_synonyms(compound) if pd.notnull(compound) else [] for compound
                                           in df["pubchem"]]
    except:
        logging.exception("No synonyms")
    df["pubchem_logp"] = [compound.xlogp if pd.notnull(compound) else np.NAN for compound in df["pubchem"]]

    return df


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
            compounds = get_pubchem_compound(smiles, "smiles")
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
    if pd.isnull(compound):
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
