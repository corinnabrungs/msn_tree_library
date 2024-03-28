import pandas as pd
import logging
from tqdm import tqdm

import pandas_utils
from date_utils import iso_datetime_now
from meta_constants import MetaColumns
from pandas_utils import (
    left_merge_retain_index,
    add_column_prefix,
    update_dataframes,
    create_missing_columns,
    read_dataframe,
)


tqdm.pandas()


def dictionary_of_np_search(df):
    from rdkit_mol_identifiers import split_inchikey

    if MetaColumns.split_inchikey not in df and MetaColumns.inchikey in df:
        df[MetaColumns.split_inchikey] = [
            split_inchikey(inchikey) for inchikey in df[MetaColumns.inchikey]
        ]
    # only merge on id column where id is notnull
    results = df[[MetaColumns.split_inchikey]][
        df[MetaColumns.split_inchikey].notnull()
    ].copy()
    if len(results) == 0:
        return df

    logging.info("Search dictionary of NP list by first block of inchikey")
    # download needed
    prefix = "dictionary_np_"
    dnp_df = read_dataframe("data/dnp_cleaned.tsv")
    if MetaColumns.split_inchikey not in dnp_df and MetaColumns.inchikey in dnp_df:
        dnp_df[MetaColumns.split_inchikey] = [
            split_inchikey(inchikey) for inchikey in dnp_df[MetaColumns.inchikey]
        ]
    # need unique split_inchikey rows for dnp to merge later
    dnp_df = dnp_df[[MetaColumns.split_inchikey]].drop_duplicates()
    dnp_df["entry"] = True

    results = left_merge_retain_index(results, dnp_df, on=MetaColumns.split_inchikey)
    results = add_column_prefix(
        results, prefix, columns_to_keep=MetaColumns.split_inchikey
    )
    results[MetaColumns.date_dictionary_np_search] = iso_datetime_now()
    return update_dataframes(results, df)
