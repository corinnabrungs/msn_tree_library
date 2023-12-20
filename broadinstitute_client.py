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


def broad_list_search(df):
    from rdkit_mol_identifiers import split_inchikey
    from drug_utils import map_clinical_phase_to_number

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

    logging.info("Search broad institute list of drugs by first block of inchikey")
    # download from: https://clue.io/repurposing#download-data
    prefix = "broad_"
    broad_df = read_dataframe("data/broad_institute_drug_list.csv")

    broad_df[MetaColumns.split_inchikey] = [
        split_inchikey(inchikey) for inchikey in broad_df["InChIKey"]
    ]
    broad_df = broad_df.drop(columns=["InChIKey"])

    # need unique split_inchikey rows for broad to merge later
    broad_df = broad_df.drop_duplicates(subset=MetaColumns.split_inchikey)

    results = left_merge_retain_index(results, broad_df, on=MetaColumns.split_inchikey)
    results = add_column_prefix(
        results, prefix, columns_to_keep=MetaColumns.split_inchikey
    )

    # converting the clinical phases to numbers (remove phase,
    # preclinic (as 0.5), or launched)
    results = create_missing_columns(results, ["broad_clinical_phase"])

    results["broad_clinical_phase"] = [
        map_clinical_phase_to_number(phase) for phase in results["broad_clinical_phase"]
    ]
    # results = pandas_utils.make_str_floor_to_int_number(
    #     df,
    #     ["broad_clinical_phase"],
    # )

    results[MetaColumns.date_broad_drug_list] = iso_datetime_now()
    return update_dataframes(results, df)
