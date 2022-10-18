import numpy as np
import pandas as pd
from datetime import date
import re

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


def ensure_synonyms_column(df: pd.DataFrame) -> pd.DataFrame:
    df["synonyms"] = df.apply(lambda row: get_all_synonyms(row), axis=1)
    return df

def get_all_synonyms(row):
    synonyms = [
        get_or_else(row, "Product Name"),
        get_or_else(row, "CAS No."),
        get_or_else(row, "CAS"),
    ]

    old = row["synonyms"] if "synonyms" in row else None
    if isinstance(old, str):
        synonyms.append(old)
    elif old is not None:
        synonyms = synonyms + old

    synonyms.extend([s.strip() for s in str(get_or_else(row, "Synonyms", "")).split(";")])

    synonyms = [x.strip() for x in synonyms if x]
    seen = set()
    unique = [x for x in synonyms if x.lower() not in seen and not seen.add(x.lower())]
    return unique

def get_or_else(row, key, default=None):
    return row[key] if key in row and not pd.isnull(row[key]) else default


def cleanup_file(metadata_file, query_pubchem: bool = True, calc_identifiers: bool = True, pubchem_search: bool = True, query_chembl: bool = True):
    logging.info("Will run on %s", metadata_file)

    # import df
    if metadata_file.endswith(".tsv"):
        df = pd.read_csv(metadata_file, sep="\t")
    else:
        df = pd.read_csv(metadata_file, sep=",")

    ensure_synonyms_column(df)

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

    # get PubChem information based on inchikey, smiles, and Inchi
    if pubchem_search:
        logging.info("Search PubChem by structure")
        df = pubchem_search_by_structure(df)

    # extract ids like the UNII, ...
    df = ensure_synonyms_column(df)
    df = extract_synonym_ids(df)

    # get ChEMBL information based on inchikey
    if query_chembl:
        logging.info("Search ChEMBL by chemblid or inchikey")
        df = chembl_search(df)


    # drop mol
    df = df.drop('mol', axis=1)

    # export metadata file
    out_file = add_suffix(metadata_file, "cleaned")
    logging.info("Exporting to file %s", out_file)
    if metadata_file.endswith(".tsv"):
        df.to_csv(out_file, sep="\t", index=False)
    else:
        df.to_csv(out_file, sep=",", index=False)


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
        dfa = df[["Cat. No.", "Product Name", "synonyms", "CAS No.", "Smiles", "pubchem_cid", "isomeric_smiles",
                  "canonical_smiles", "lib_plate_well", "URL", "Target", "Information", "Pathway", "Research Area",
                  "Clinical Information"]].copy()

        dfa["Source"] = "MCE"
        dfb = df[["Cat. No.", "Product Name", "synonyms", "CAS No.", "pubchem_cid", "isomeric_smiles",
                  "canonical_smiles", "lib_plate_well", "URL", "Target", "Information", "Pathway", "Research Area",
                  "Clinical Information"]].copy()
        dfb["Smiles"] = dfb["isomeric_smiles"]
        dfb["Source"] = "PubChem"

        df = pd.concat([dfb, dfa], ignore_index=True, sort=False)

    return df[df["Smiles"].notna()]

def pubchem_search_by_structure(df) -> pd.DataFrame:
    compounds = [client.search_pubchem_by_structure(smiles, inchi, inchikey) for inchikey, smiles, inchi in zip(df["inchi_key"], df["Smiles"], df["inchi"])]

    df["pubchem_cid_parent"] = pd.array([compound.cid if not pd.isnull(compound) else np.NAN for compound in compounds],
                                dtype=pd.Int64Dtype())
    df["compound_name"] = [compound.synonyms[0] if not pd.isnull(compound) else np.NAN for compound in compounds]
    df["iupac"] = [compound.iupac_name if not pd.isnull(compound) else np.NAN for compound in compounds]
    df["synonyms"] = df["synonyms"] + [compound.synonyms if not pd.isnull(compound) else [] for compound in compounds]
    df["pubchem_logp"] = [compound.xlogp if not pd.isnull(compound) else np.NAN for compound in compounds]

    return df

def chembl_search(df) -> pd.DataFrame:
    compounds = [client.get_chembl_mol(chembl_id, inchi_key) for chembl_id, inchi_key in zip(df["chembl_id"], df["inchi_key"])]

    df["chembl_id"] = [compound["molecule_chembl_id"] if not pd.isnull(compound) else np.NAN for compound in compounds]
    # df["compound_name"] = df["compound_name"] + [compound["pref_name"] if not pd.isnull(compound) else np.NAN for compound in compounds]
    df["molecular_species"] = [compound["molecule_properties"]["molecular_species"] if not pd.isnull(compound) else np.NAN for compound in compounds]
    df["prodrug"] = [compound["prodrug"] if not pd.isnull(compound) else np.NAN for compound in compounds]
    df["clinical_phase"] = [compound["max_phase"] if not pd.isnull(compound) else np.NAN for compound in compounds]
    df["first_approval"] = pd.array([compound["first_approval"] if not pd.isnull(compound) else np.NAN for compound in compounds], dtype=pd.Int64Dtype())
    df["withdrawn"] = [compound["withdrawn_flag"] if not pd.isnull(compound) else np.NAN for compound in compounds]
    df["withdrawn_class"] = [compound["withdrawn_class"] if not pd.isnull(compound) else np.NAN for compound in compounds]
    df["withdrawn_reason"] = [compound["withdrawn_reason"] if not pd.isnull(compound) else np.NAN for compound in compounds]
    df["withdrawn_year"] = pd.array([compound["withdrawn_year"] if not pd.isnull(compound) else np.NAN for compound in compounds], dtype=pd.Int64Dtype())
    df["withdrawn_country"] = [compound["withdrawn_country"] if not pd.isnull(compound) else np.NAN for compound in compounds]
    df["oral"] = [compound["oral"] if not pd.isnull(compound) else np.NAN for compound in compounds]
    df["parenteral"] = [compound["parenteral"] if not pd.isnull(compound) else np.NAN for compound in compounds]
    df["topical"] = [compound["topical"] if not pd.isnull(compound) else np.NAN for compound in compounds]
    df["natural_product"] = [compound["natural_product"] if not pd.isnull(compound) else np.NAN for compound in compounds]
    df["usan_stem_definition"] = [compound["usan_stem_definition"] if not pd.isnull(compound) else np.NAN for compound in compounds]
    df["chembl_alogp"] = [compound["molecule_properties"]["alogp"] if not pd.isnull(compound) else np.NAN for compound in compounds]
    df["chembl_clogp"] = [compound["molecule_properties"]["cx_logp"] if not pd.isnull(compound) else np.NAN for compound in compounds]

    ## dont overwrite
    # df["synonyms"] = df["synonyms"] + [compound["molecule_synonyms"] if not pd.isnull(compound) else [] for compound in compounds]
    df["indication"] = [compound["indication_class"] if not pd.isnull(compound) else np.NAN for compound in compounds]

    return df


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
    drugbank_generator = (name.upper() for name in synonyms if name.upper().startswith("DB-"))
    return next(drugbank_generator, None)

def extract_synonym_ids(df: pd.DataFrame) -> pd.DataFrame:
    df["unii"] = [find_unii(synonyms) for synonyms in df["synonyms"]]
    df["schembl_id"] = [find_schembl(synonyms) for synonyms in df["synonyms"]]
    df["chembl_id"] = [find_chembl_id(synonyms) for synonyms in df["synonyms"]]
    df["zinc_id"] = [find_zinc(synonyms) for synonyms in df["synonyms"]]
    df["drugbank_id"] = [find_drugbank(synonyms) for synonyms in df["synonyms"]]

    return df


if __name__ == "__main__":
    cleanup_file(r"data\test_metadata.tsv", query_pubchem=True)
    # cleanup_file("data\lib_formatted_pubchem_mce.tsv", query_pubchem=True)
    # cleanup_file("data\mce_library_add_compounds.tsv", query_pubchem=True)
