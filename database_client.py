import json

import pandas as pd
import numpy as np
from pubchempy import get_compounds, Compound
from chembl_webresource_client.new_client import new_client as chembl
import requests
import logging

# openfda by name (uppercase)
OPENFDA_URL = r"https://api.fda.gov/other/substance.json?search=names.name:%22{}%22"

#opfenfda by unii
    # r"https://api.fda.gov/other/substance.json?search=unii:"{}""
#openfda by substance name
    # r"https://api.fda.gov/drug/drugsfda.json?search=openfda.substance_name:%22{}%22&limit=1"
#openfda by cas
    # r"https://api.fda.gov/other/substance.json?search=codes.code:"{}""


logging.getLogger('pubchempy').setLevel(logging.DEBUG)


def pubchem_compound_score(comp: Compound):
    smiles = comp.canonical_smiles
    if not smiles:
        return 0
    return 1000 - str(smiles).count(".")


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
        raise ValueError("At least one structure identifier needs to be a value")

    compounds = None

    if inchikey:
        compounds = get_compounds(inchikey, "inchikey")
    if not compounds and smiles:
        compounds = get_compounds(smiles, "smiles")
    if not compounds and inchi:
        compounds = get_compounds(inchi, "inchi")
    if not compounds:
        logging.info("NO PUBCHEM FOR: smiles:{}  inchi:{}   inchikey:{}".format(smiles, inchi, inchikey))
        return None
    else:
        compounds.sort(key=lambda comp: pubchem_compound_score(comp), reverse=True)
        return compounds[0]


def get_chembl_mol_by_inchikey(inchi_key):
    try:
        compounds = chembl.molecule.filter(molecule_structures__standard_inchi_key=inchi_key)
        if not compounds:
            logging.info("NO ChEMBL FOR: inchikey:{}".format(inchi_key))
            return None
        else:
            return compounds[0]
    except Exception as e:
        logging.warning("Error during chembl query:",e)
        return None


def get_openfda_information(name):
    url = OPENFDA_URL.format(name)
    response = requests.get(url)
    response.raise_for_status()
    return response.json()

# def get drugcentral_information(name):
