import numpy as np
import pandas as pd
from datetime import date
import re
from tqdm import tqdm

from molmass import Formula
from rdkit.Chem import Descriptors
from rdkit.Chem import AllChem as Chem
from rdkit.Chem.MolStandardize import rdMolStandardize
from pathlib import Path
import logging
import mol_identifiers as molid
import database_client as client
import drugcentral_postgresql_query as drugcentral_query

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

# query braod list
def map_clinical_phase_to_number(phase):
    match (str(phase)):
        case "" | "None" | "NaN" | "nan" | "np.nan" | "0" | "No Development Reported":
            return 0
        case "Preclinical":
            return 0.5
        case "1" | "1.0" | "Phase 1":
            return 1
        case "2" | "2.0" | "Phase 2":
            return 2
        case "3" | "3.0" | "Phase 3":
            return 3
        case "4" | "4.0" | "Phase 4" | "Launched":
            return 4
        case _:
            return phase

def get_clinical_phase_description(number):
    match (str(number)):
        case "" | "None" | "NaN" | "nan" | "np.nan" | "0" | "0.0":
            return ""
        case "0.5":
            return "Preclinical"
        case "1.0" | "1":
            return "Phase 1"
        case "2.0" | "2":
            return "Phase 2"
        case "3.0" | "3":
            return "Phase 3"
        case "4.0" | "4":
            return "Launched"
        case _:
            return number

def broad_list_search(df):
    # download from: https://clue.io/repurposing#download-data
    prefix = "broad_"
    broad_df = pd.read_csv("data/broad_institute_drug_list.csv")
    broad_df = broad_df.add_prefix(prefix)

    if "split_inchi_key" not in df and "inchi_key" in df:
        df["split_inchi_key"] = [str(inchikey).split("-")[0] if inchikey is not None else None for inchikey in df['inchi_key']]
    broad_df["split_inchi_key"] = [str(inchikey).split("-")[0] for inchikey in broad_df["{}InChIKey".format(prefix)]]

    merged_df = pd.merge(df, broad_df, on="split_inchi_key", how="left")
    # converting the clinical phases (from broad institute, chembl, provider, or else) to numbers (remove phase, preclinic (as 0.5), or launched)
    merged_df["broad_clinical_phase"] = [map_clinical_phase_to_number(phase) for phase in
                                         merged_df["broad_clinical_phase"]]
    merged_df["clinical_phase"] = [map_clinical_phase_to_number(phase) for phase in merged_df["clinical_phase"]]
    merged_df["Clinical Information"] = [map_clinical_phase_to_number(phase) for phase in
                                         merged_df["Clinical Information"]]
    # Comparing the clinical phases, only store the highest number in clinical phase column
    merged_df["clinical_phase"] = merged_df[['broad_clinical_phase', 'clinical_phase', 'Clinical Information']].max(
        axis=1)
    # Converting numbers back to phase X, launched or preclinic
    merged_df["clinical_phase_description"] = [get_clinical_phase_description(number) for number in
                                               merged_df["clinical_phase"]]
    return merged_df.drop(["broad_clinical_phase", "Clinical Information", "{}InChIKey".format(prefix)], axis=1)


def find_in_drugbank(drugbank_df, row):
    if pd.notnull(row["drugbank_id"]):
        return row["drugbank_id"]

    dbid = None
    # pubchem id first, then CHEMBL, then synonyms
    if row["inchi_key"]:
        dbid = next((d for d in drugbank_df[drugbank_df["inchi_key"]==row["inchi_key"]]["drugbank_id"]), None)
    if not dbid and row["pubchem_cid_parent"]:
        dbid = next((d for d in drugbank_df[drugbank_df["pubchem_cid"]==row["pubchem_cid_parent"]]["drugbank_id"]), None)
    if not dbid and row["chembl_id"]:
        dbid = next((d for d in drugbank_df[drugbank_df["chembl_id"]==row["chembl_id"]]["drugbank_id"]), None)
    if not dbid and row["unii"]:
        dbid = next((d for d in drugbank_df[drugbank_df["unii"]==row["unii"]]["drugbank_id"]), None)
    if not dbid and row["CAS No."]:
        dbid = next((d for d in drugbank_df[drugbank_df["cas"]==row["CAS No."]]["drugbank_id"]), None)
    if not dbid and row["split_inchi_key"]:
        dbid = next((d for d in drugbank_df[drugbank_df["split_inchi_key"]==row["split_inchi_key"]]["drugbank_id"]), None)
    # if not dbid and row["synonyms"]:
    #     dbid = next((d for d in drugbank_df[drugbank_df["name"] in row["synonyms"]]["drugbank_id"]), None)
    return dbid

def drugbank_list_search(df):
    # download from: https://go.drugbank.com/releases/latest, approved access needed, xml extraction to tsv by drugbank_extraction.py
    prefix = "drugbank_"
    drugbank_df = pd.read_csv("data/drugbank.tsv", sep="\t")
    drugbank_df["pubchem_cid"] = pd.array(drugbank_df["pubchem_cid"], dtype=pd.Int64Dtype())

    if "split_inchi_key" not in df and "inchi_key" in df:
        df["split_inchi_key"] = [str(inchikey).split("-")[0] if inchikey is not None else None for inchikey in df['inchi_key']]
    drugbank_df["split_inchi_key"] = [str(inchikey).split("-")[0] for inchikey in drugbank_df["inchi_key"]]

    df["drugbank_id"] = None
    # find drugbank IDs in drugbank table by PubChem, ChEMBL etc
    df["drugbank_id"] = df.apply(lambda row: find_in_drugbank(drugbank_df, row), axis=1)

    drugbank_df = drugbank_df.add_prefix(prefix)
    merged_df = pd.merge(df, drugbank_df, left_on="drugbank_id", right_on="drugbank_drugbank_id", how="left")
    return merged_df.drop(["drugbank_drugbank_id", "{}inchi_key".format(prefix), "{}smiles".format(prefix), "{}split_inchi_key".format(prefix)], axis=1)

def drugcentral_search(df):
    if "split_inchi_key" not in df and "inchi_key" in df:
        df["split_inchi_key"] = [str(inchikey).split("-")[0] if pd.notnull(inchikey) else None for inchikey in df['inchi_key']]
    try:
        drugcentral_query.connect()
        logging.info("Searching in DrugCentral")
        results = [drugcentral_query.drugcentral_postgresql(inchikey, split_inchikey) for inchikey, split_inchikey in
                 tqdm(zip(df["inchi_key"], df["split_inchi_key"]))]
        logging.info("DrugCentral search done")

        # row[1] is the data row[0] is the columns
        first_entry_columns = next((row[0] for row in results), [])
        columns = [col.name for col in first_entry_columns]
        elements = len(columns)
        data = [row[1] if pd.notnull(row[1]) else (None, )*elements for row in results]
        dc_df = pd.DataFrame(data=data, columns=columns, index=df.index)

        return pd.concat([df, dc_df], axis=1)

    finally:
        drugcentral_query.deconnect()


def cleanup_file(metadata_file, query_pubchem: bool = True, calc_identifiers: bool = True, pubchem_search: bool = True,
                 query_chembl: bool = True, query_broad_list=False, query_drugbank_list=False, query_drugcentral=False):
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

    # get broad institute information based on split inchikey
    if query_broad_list:
        logging.info("Search broad institute list of drugs by first block of inchikey")
        df = broad_list_search(df)

        # get broad institute information based on split inchikey
    if query_drugbank_list:
        logging.info("Search drugbank list by inchikey, pubchem_id, chembl_id, cas, split inchikey, etc.")
        df = drugbank_list_search(df)

    if query_drugcentral:
        logging.info("Search drugcentral by inchikey")
        df = drugcentral_search(df)

    # drop mol
    df = df.drop('mol', axis=1)

    # export metadata file
    out_file = add_suffix(metadata_file, "cleaned")
    logging.info("Exporting to file %s", metadata_file)
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
    df["split_inchi_key"] = [str(inchikey).split("-")[0] for inchikey in df['inchi_key']]
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


if __name__ == "__main__":
    cleanup_file(r"data\test_metadata_small.tsv", query_pubchem=True, query_broad_list=True, query_drugbank_list=True,
                 query_drugcentral=True)
    # cleanup_file("data\lib_formatted_pubchem_mce.tsv", query_pubchem=True)
    # cleanup_file("data\mce_library_add_compounds.tsv", query_pubchem=True)
