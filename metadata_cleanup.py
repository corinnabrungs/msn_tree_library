import numpy as np
import pandas as pd
from datetime import date

from molmass import Formula
from rdkit.Chem import Descriptors
from rdkit.Chem import AllChem as Chem
from rdkit.Chem.MolStandardize import rdMolStandardize
from pathlib import Path
import logging
import mol_identifiers as molid
import database_client as client

logging.basicConfig(format="%(asctime)s - %(message)s", level=logging.DEBUG)


def add_suffix(filename, suffix):
    p = Path(filename)
    return "{0}_{1}{2}".format(Path.joinpath(p.parent, p.stem), suffix, p.suffix)


def cleanup_file(metadata_file, query_pubchem: bool = True, calc_identifiers: bool = True, pubchem_search: bool = True):
    logging.info("Will run on %s", metadata_file)

    # import df
    if metadata_file.endswith(".tsv"):
        df = pd.read_csv(metadata_file, sep="\t")
    else:
        df = pd.read_csv(metadata_file, sep=",")

    # Query pubchem by name and CAS
    if query_pubchem:
        logging.info("Search PubChem by name")
        df = pubchem_search_structure_by_name(df)

    # get mol from smiles or inchi
    # calculate all identifiers from mol - exact_mass, ...
    if calc_identifiers:
        add_molid_columns(df)
    # drop duplicates
    df = df.drop_duplicates(['Product Name', 'lib_plate_well', "inchi_key"], keep="first").sort_index()

    # get PubChem information based on canonical_smiles and Inchi
    if pubchem_search:
        logging.info("Search PubChem by structure")
        pubchem_search_by_structure(df)



    # drop mol
    df = df.drop('mol', axis=1)

    # export metadata file
    out_file = add_suffix(metadata_file, "cleaned")
    logging.info("Exporting to file %s", out_file)
    if metadata_file.endswith(".tsv"):
        df.to_csv(out_file, sep="\t")
    else:
        df.to_csv(out_file, sep=",")


def add_molid_columns(df):
    # first strip any salts
    df["Smiles"] = [molid.split_smiles_major_mol(smiles) if not pd.isnull(smiles) else np.NAN for smiles in df["Smiles"]]
    df["mol"] = [Chem.MolFromSmiles(smiles) if not pd.isnull(smiles) else np.NAN for smiles in df["Smiles"]]
    df["mol"] = [molid.chembl_standardize_mol(mol) if not pd.isnull(mol) else np.NAN for mol in df["mol"]]
    df["canonical_smiles"] = [molid.mol_to_canon_smiles(mol) for mol in df["mol"]]
    df["Smiles"] = [molid.mol_to_isomeric_smiles(mol) for mol in df["mol"]]
    df["exact_mass"] = [molid.exact_mass_from_mol(mol) for mol in df["mol"]]
    df["inchi"] = [molid.inchi_from_mol(mol) for mol in df["mol"]]
    df["inchi_key"] = [molid.inchikey_from_mol(mol) for mol in df["mol"]]
    df["formula"] = [molid.formula_from_mol(mol) for mol in df["mol"]]
    return df


def pubchem_search_structure_by_name(df) -> pd.DataFrame:
    compounds = [client.search_pubchem_by_name(str(cas)) if not pd.isnull(cas) else np.NAN for cas in df["CAS No."]]
    compounds = [client.search_pubchem_by_name(str(name)) if pd.isnull(comp) else comp for comp, name in
                 zip(compounds, df["Product Name"])]
    # only one compound was found as CAS-
    compounds = [client.search_pubchem_by_name("CAS-{}".format(cas)) if pd.isnull(comp) else comp for comp, cas in
                 zip(compounds, df["CAS No."])]

    df["pubchem_cid"] = pd.array([compound.cid if not pd.isnull(compound) else np.NAN for compound in compounds],
                                dtype=pd.Int64Dtype())
    df["isomeric_smiles"] = [compound.isomeric_smiles if not pd.isnull(compound) else np.NAN for compound in compounds]
    df["canonical_smiles"] = [compound.canonical_smiles if not pd.isnull(compound) else np.NAN for compound in
                              compounds]

    # concat the new structures with the old ones
    if "Smiles" in df:
        dfa = df[["Cat. No.", "Product Name", "Synonyms", "CAS No.", "Smiles", "pubchem_cid", "isomeric_smiles",
                  "canonical_smiles", "lib_plate_well", "URL", "Target", "Information", "Pathway", "Research Area",
                  "Clinical Information"]].copy()

        dfa["Source"] = "MCE"
        dfb = df[["Cat. No.", "Product Name", "Synonyms", "CAS No.", "pubchem_cid", "isomeric_smiles",
                  "canonical_smiles", "lib_plate_well", "URL", "Target", "Information", "Pathway", "Research Area",
                  "Clinical Information"]].copy()
        dfb["Smiles"] = dfb["isomeric_smiles"]
        dfb["Source"] = "PubChem"

        df = pd.concat([dfb, dfa], ignore_index=True, sort=False)

    return df[df["Smiles"].notna()]

def pubchem_search_by_structure(df) -> pd.DataFrame:
    compounds = [client.search_pubchem_by_structure(inchikey, smiles, inchi) for inchikey, smiles, inchi in zip(df["inchi_key"], df["Smiles"], df["inchi"])]

    df["pubchem_cid_parent"] = pd.array([compound.cid if not pd.isnull(compound) else np.NAN for compound in compounds],
                                dtype=pd.Int64Dtype())
    df["iupac"] = [compound.iupac_name if not pd.isnull(compound) else np.NAN for compound in compounds]


if __name__ == "__main__":
    cleanup_file(r"data\test_metadata.tsv", query_pubchem=True)
    # cleanup_file("data\lib_formatted_pubchem_mce.tsv", query_pubchem=True)
    # cleanup_file("data\mce_library_add_compounds.tsv", query_pubchem=True)
