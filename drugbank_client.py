import logging

import pandas as pd

import pandas_utils
from date_utils import iso_datetime_now
from meta_constants import MetaColumns
from rdkit_mol_identifiers import split_inchikey
from pandas_utils import (
    isnull,
    notnull,
    create_missing_columns,
    combine_dfs_fill_missing_values,
    make_str_floor_to_int_number,
    update_dataframes,
    add_column_prefix,
    left_merge_retain_index,
    notnull_not_empty,
)
from drug_utils import map_clinical_phase_to_number
from tqdm import tqdm

tqdm.pandas()


def map_drugbank_approval(status):
    match (str(status)):
        case "approved" | "withdrawn":
            return 4
        case _:
            return None


def find_in_drugbank(drugbank_df, row):
    if notnull_not_empty(row["drugbank_id"]):
        return row["drugbank_id"]

    dbid = None
    # pubchem id first, then CHEMBL, then synonyms
    if notnull_not_empty(row["drugbank_id"]):
        dbid = next(
            (
                d
                for d in drugbank_df[drugbank_df["drugbank_id"] == row["drugbank_id"]][
                    "drugbank_id"
                ]
            ),
            None,
        )
    if isnull(dbid) and notnull_not_empty(row["inchikey"]):
        dbid = next(
            (
                d
                for d in drugbank_df[drugbank_df["inchikey"] == row["inchikey"]][
                    "drugbank_id"
                ]
            ),
            None,
        )
    if isnull(dbid) and notnull_not_empty(row[MetaColumns.pubchem_cid]):
        dbid = next(
            (
                d
                for d in drugbank_df[
                    drugbank_df["pubchem_cid"] == row[MetaColumns.pubchem_cid]
                ]["drugbank_id"]
            ),
            None,
        )
    if isnull(dbid) and notnull_not_empty(row[MetaColumns.input_pubchem_cid]):
        dbid = next(
            (
                d
                for d in drugbank_df[
                    drugbank_df["pubchem_cid"] == row[MetaColumns.input_pubchem_cid]
                ]["drugbank_id"]
            ),
            None,
        )
    if isnull(dbid) and notnull_not_empty(row["chembl_id"]):
        dbid = next(
            (
                d
                for d in drugbank_df[drugbank_df["chembl_id"] == row["chembl_id"]][
                    "drugbank_id"
                ]
            ),
            None,
        )
    if isnull(dbid) and notnull_not_empty(row["unii"]):
        dbid = next(
            (d for d in drugbank_df[drugbank_df["unii"] == row["unii"]]["drugbank_id"]),
            None,
        )
    if isnull(dbid) and "cas" in row and notnull(row["cas"]):
        dbid = next(
            (d for d in drugbank_df[drugbank_df["cas"] == row["cas"]]["drugbank_id"]),
            None,
        )
    if isnull(dbid) and notnull_not_empty(row["split_inchikey"]):
        dbid = next(
            (
                d
                for d in drugbank_df[
                    drugbank_df["split_inchikey"] == row["split_inchikey"]
                ]["drugbank_id"]
            ),
            None,
        )
    if isnull(dbid) and notnull_not_empty(row["synonyms"]):
        lower_synonyms = [s.lower().strip() for s in row["synonyms"] if notnull(s)]
        filtered_dbank_df = drugbank_df[
            drugbank_df["lower_case_compound_name"].isin(lower_synonyms)
        ]
        dbid = next((d for d in filtered_dbank_df["drugbank_id"]), None)
    return dbid


def drugbank_list_search(df: pd.DataFrame):
    # download from: https://go.drugbank.com/releases/latest, approved access needed, xml extraction to tsv by
    # drugbank_extraction.py
    prefix = "drugbank_"
    drugbank_df = pd.read_csv("data/drugbank.tsv", sep="\t")
    drugbank_df["lower_case_compound_name"] = [
        cname.lower().strip() if notnull(cname) else None
        for cname in drugbank_df["compound_name"]
    ]
    drugbank_df = make_str_floor_to_int_number(drugbank_df, MetaColumns.pubchem_cid)

    if "split_inchikey" not in df and "inchikey" in df:
        df["split_inchikey"] = [split_inchikey(inchikey) for inchikey in df["inchikey"]]
    drugbank_df["split_inchikey"] = [
        split_inchikey(inchikey) for inchikey in drugbank_df["inchikey"]
    ]

    # find drugbank IDs in drugbank table by PubChem, ChEMBL etc
    df["drugbank_id"] = df.progress_apply(
        lambda row: find_in_drugbank(drugbank_df, row), axis=1
    )

    # TODO check issue that sometimes we get series not hashable exception
    df["drugbank_id"] = [
        did if isinstance(did, str) else None for did in df["drugbank_id"]
    ]

    # only merge on id column where id is notnull
    results = df[[MetaColumns.drugbank_id]][
        df[MetaColumns.drugbank_id].notnull()
    ].copy()
    if len(results) == 0:
        return df

    results = left_merge_retain_index(results, drugbank_df, on=MetaColumns.drugbank_id)
    results = results.drop(
        columns=[
            "inchikey",
            "smiles",
            "inchi",
            MetaColumns.isomeric_smiles,
            MetaColumns.canonical_smiles,
            "split_inchikey",
        ],
        errors="ignore",
    )

    # rename only a few columns with prefix
    results = add_column_prefix(
        results,
        prefix,
        columns_to_keep=[
            "unii",
            "chembl_id",
            "pubchem_cid",
            "compound_name",
            "cas",
            MetaColumns.drugbank_id,
        ],
    )

    # fill NA values in original data with drugbank data
    return update_dataframes(results, df)


def drugbank_search_add_columns(df):
    logging.info(
        "Search drugbank list by drugbank_id, inchikey, pubchem_id, chembl_id, cas, split inchikey, etc."
    )

    df = drugbank_list_search(df)
    df = create_missing_columns(df, ["drugbank_approved"])
    df["drugbank_clinical_phase"] = [
        map_drugbank_approval(status) for status in df["drugbank_approved"]
    ]
    df["drugbank_clinical_phase"] = [
        map_clinical_phase_to_number(phase) for phase in df["drugbank_clinical_phase"]
    ]

    df[MetaColumns.date_drugbank_search] = iso_datetime_now()
    return df
