import logging

import pandas as pd

import pandas_utils
from rdkit_mol_identifiers import split_inchikey
from pandas_utils import isnull, notnull
from tqdm import tqdm

tqdm.pandas()


def map_drugbank_approval(status):
    match (str(status)):
        case "approved" | "withdrawn":
            return 4
        case _:
            return None


def find_in_drugbank(drugbank_df, row):
    if notnull(row["drugbank_id"]):
        return row["drugbank_id"]

    dbid = None
    # pubchem id first, then CHEMBL, then synonyms
    if notnull(row["drugbank_id"]):
        dbid = next((d for d in drugbank_df[drugbank_df["drugbank_id"] == row["drugbank_id"]]["drugbank_id"]), None)
    if notnull(row["inchikey"]):
        dbid = next((d for d in drugbank_df[drugbank_df["inchikey"] == row["inchikey"]]["drugbank_id"]), None)
    if isnull(dbid) and not isnull(row["pubchem_cid_parent"]):
        dbid = next((d for d in drugbank_df[drugbank_df["pubchem_cid"] == row["pubchem_cid_parent"]]["drugbank_id"]),
                    None)
    if isnull(dbid) and notnull(row["chembl_id"]):
        dbid = next((d for d in drugbank_df[drugbank_df["chembl_id"] == row["chembl_id"]]["drugbank_id"]), None)
    if isnull(dbid) and notnull(row["unii"]):
        dbid = next((d for d in drugbank_df[drugbank_df["unii"] == row["unii"]]["drugbank_id"]), None)
    if isnull(dbid) and "cas" in row and notnull(row["cas"]):
        dbid = next((d for d in drugbank_df[drugbank_df["cas"] == row["cas"]]["drugbank_id"]), None)
    if isnull(dbid) and notnull(row["split_inchikey"]):
        dbid = next((d for d in drugbank_df[drugbank_df["split_inchikey"] == row["split_inchikey"]]["drugbank_id"]),
                    None)
    # if isnull(dbid) and row["synonyms"]:
    #     dbid = next((d for d in drugbank_df[drugbank_df["compound_name"] in row["synonyms"]]["drugbank_id"]), None)
    return dbid


def drugbank_list_search(df):
    # download from: https://go.drugbank.com/releases/latest, approved access needed, xml extraction to tsv by
    # drugbank_extraction.py
    prefix = "drugbank_"
    drugbank_df = pd.read_csv("data/drugbank.tsv", sep="\t")
    # TODO check and maybe use str
    drugbank_df["pubchem_cid"] = pd.array(drugbank_df["pubchem_cid"], dtype=pd.Int64Dtype())

    if "split_inchikey" not in df and "inchikey" in df:
        df["split_inchikey"] = [split_inchikey(inchikey) for inchikey in df['inchikey']]
    drugbank_df["split_inchikey"] = [split_inchikey(inchikey) for inchikey in drugbank_df["inchikey"]]

    # find drugbank IDs in drugbank table by PubChem, ChEMBL etc
    df["drugbank_id"] = df.progress_apply(lambda row: find_in_drugbank(drugbank_df, row), axis=1)

    # TODO check issue that sometimes we get series not hashable exception
    df["drugbank_id"] = [did if isinstance(did, str) else None for did in df["drugbank_id"]]

    drugbank_df = drugbank_df.add_prefix(prefix)
    merged_df = pd.merge(df, drugbank_df, left_on="drugbank_id", right_on="drugbank_drugbank_id", how="left")

    filldf = merged_df[
        ["drugbank_unii", "drugbank_chembl_id", "drugbank_pubchem_cid", "drugbank_compound_name", "drugbank_inchikey",
         "drugbank_split_inchikey", "drugbank_cas"]].copy()
    filldf.columns = filldf.columns.str.removeprefix(prefix)
    merged_df = pandas_utils.combine_dfs_fill_missing_values(merged_df, filldf)

    # merged_df["unii"].combine_first(merged_df["{}unii".format(prefix)])
    # merged_df["chembl_id"].combine_first(merged_df["{}chembl_id".format(prefix)], inplace=True)
    # merged_df.combine_first(merged_df["{}name".format(prefix)], inplace=True)
    return merged_df.drop(
        columns=["drugbank_drugbank_id", "{}inchikey".format(prefix), "{}smiles".format(prefix), f"{prefix}unii",
                 f"{prefix}chembl_id", "{}split_inchikey".format(prefix)], errors="ignore")


def drugbank_search_add_columns(df):
    logging.info("Search drugbank list by drugbank_id, inchikey, pubchem_id, chembl_id, cas, split inchikey, etc.")
    df = drugbank_list_search(df)
    if "drugbank_approved" in df.columns:
        df["drugbank_approved_number"] = [map_drugbank_approval(status) for status in df["drugbank_approved"]]
    else:
        df["drugbank_approved_number"] = None

    df["clinical_phase"] = df[df.columns[df.columns.isin(["clinical_phase", "drugbank_approved_number"])]].max(axis=1)
    df["any_phase"] = df["drugbank_approved"].notna() | (df["clinical_phase"] > 0)
    return df
