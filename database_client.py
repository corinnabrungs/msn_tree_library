import pandas as pd
import numpy as np
from pubchempy import get_compounds, Compound
from chembl_webresource_client.new_client import new_client as chembl
import logging
logging.getLogger('pubchempy').setLevel(logging.DEBUG)


def compound_score(comp: Compound):
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
        compounds.sort(key=lambda comp: compound_score(comp), reverse=True)
        return compounds[0]


def search_pubchem_by_structure(smiles, inchi, inchikey) -> Compound | None:
    """
    In pubchem many entries contain the cas as an alternative name - so searching for cas in name works often

    :param smiles:
    :param inchi:
    :param inchikey:
    :return: first compound or None
    """
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
        compounds.sort(key=lambda comp: compound_score(comp), reverse=True)
        return compounds[0]


def get_chembl_mol_by_inchikey(inchi_key):
    try:
        return chembl.molecule.filter(molecule_structures__standard_inchi_key=inchi_key)
    except Exception as e:
        logging.warning("Error during chembl query:",e)
        return None