import json

import pandas as pd
import numpy as np
from pubchempy import get_compounds, Compound
from chembl_webresource_client.new_client import new_client as chembl
import requests
import logging
# from diskcache import Cache
# cache = Cache("tmpcache")
from joblib import Memory
memory = Memory("memcache")


# openfda by name (uppercase)
OPENFDA_URL = r"https://api.fda.gov/other/substance.json?search=names.name:%22{}%22"

# opfenfda by unii
OPENFDA_UNII_URL = r"https://api.fda.gov/other/substance.json?search=unii:{}"
# openfda by substance name
# r"https://api.fda.gov/drug/drugsfda.json?search=openfda.substance_name:%22{}%22&limit=1"
# openfda by cas
# r"https://api.fda.gov/other/substance.json?search=codes.code:"{}""


# drugcentral
DRUGCENTRAL_URL = r"https://pharos-api.newdrugtargets.org/drugcentral?name={}%20"

logging.getLogger('pubchempy').setLevel(logging.DEBUG)


def pubchem_compound_score(comp: Compound):
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
        logging.info("cas:{} had NO entries".format(name_or_cas))
        return None
    else:
        compounds.sort(key=lambda comp: pubchem_compound_score(comp), reverse=True)
        return compounds[0]


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
    try:
        return get_compounds(value, key)
    except:
        logging.warning("FAILED PUBCHEM FOR: {} (as {})".format(value, key))
        return None


@memory.cache
def pubchem_get_synonyms(compound, try_n=1, max_tries=10):
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
        if try_n<max_tries:
            pubchem_get_synonyms(compound, try_n=try_n+1, max_tries=max_tries)
        else:
            logging.exception("Failed to retrieve synonyms for compound {}".format(compound.cid))
            return []



def get_chembl_mol(chembl_id=None, inchi_key=None):
    try:
        if not (chembl_id or inchi_key):
            raise ValueError("At least one chembl identifier needs to be a value")


        if chembl_id:
            comp = chembl.molecule.get(chembl_id)
            if comp:
                return comp

        compounds = None
        if not compounds and inchi_key:
            compounds = chembl.molecule.filter(molecule_structures__standard_inchi_key=inchi_key)
        if not compounds:
            logging.info("NO ChEMBL FOR: chemblid: {} or inchikey: {}".format(chembl_id, inchi_key))
            return None
        else:
            return compounds[0]
    except Exception as e:
        logging.warning("Error during chembl query:", e)
        return None

def get_openfda_information(name):
    url = OPENFDA_URL.format(name)
    response = requests.get(url)
    response.raise_for_status()
    return response.json()

def get_openfda_unii_information(unii):
    url = OPENFDA_UNII_URL.format(unii)
    response = requests.get(url)
    response.raise_for_status()
    return response.json()

def get_drugcentral_information(name):
    url = DRUGCENTRAL_URL.format(name)
    response = requests.get(url)
    response.raise_for_status()
    return response.json()

