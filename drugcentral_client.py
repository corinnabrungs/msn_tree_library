import logging
import pandas as pd
from tqdm import tqdm

import pandas_utils
import synonyms
from date_utils import iso_datetime_now
from drugcentral_postgresql_query import DRUGCENTRAL_ADD_SQL
from meta_constants import MetaColumns
from pandas_utils import notnull, update_dataframes, make_str_floor_to_int_number
from drug_utils import map_clinical_phase_to_number

import drugcentral_postgresql_query as drugcentral_query
from rdkit_mol_identifiers import split_inchikey

tqdm.pandas()


def drugcentral_search(df):
    logging.info("Search drugcentral by external identifier or inchikey")
    prefix = "drugcentral_"
    if "split_inchikey" not in df and "inchikey" in df:
        df["split_inchikey"] = [split_inchikey(inchikey) for inchikey in df["inchikey"]]
    try:
        drugcentral_query.connect()
        logging.info("Searching in DrugCentral")
        results = df.progress_apply(
            lambda row: drugcentral_query.drugcentral_for_row(row), axis=1
        )
        logging.info("DrugCentral main search done")

        # results may be all empty
        if all(v[0] is None for v in results):
            df[MetaColumns.date_drugcentral_search] = iso_datetime_now()
            return df

        dc_df = sql_results_to_df(df, results)
        # rename some columns
        dc_df = dc_df.rename(
            columns={
                "cas_reg_no": "cas",
                "name": "compound_name",
            }
        )
        # run additional query with id
        for i, sql_query in enumerate(DRUGCENTRAL_ADD_SQL):
            logging.info("Additional drugcentral query {}".format(i + 1))
            results = dc_df["id"].progress_apply(
                lambda dc_id: drugcentral_query.drugcentral_additional_query(
                    dc_id, sql_query
                )
            )
            logging.info("DrugCentral additional search done")
            add_df = sql_results_to_df(df, results)
            dc_df = update_dataframes(dc_df, add_df)

        dc_df = dc_df.add_prefix(prefix)

        dc_df[MetaColumns.date_drugcentral_search] = iso_datetime_now()
        df = update_dataframes(dc_df, df).copy()

        syncol = f"{prefix}synonyms"
        df = synonyms.add_synonyms_columns(
            df, new_syn_column_header=syncol, prepend=False
        )

        if "drugcentral_administration" in df.columns:
            df["drugcentral_clinical_phase"] = [
                4 if notnull(status) else None
                for status in df["drugcentral_administration"]
            ]
        else:
            df["drugcentral_clinical_phase"] = None

        df["drugcentral_clinical_phase"] = [
            map_clinical_phase_to_number(phase)
            for phase in df["drugcentral_clinical_phase"]
        ]
        return df
    finally:
        drugcentral_query.deconnect()


def sql_results_to_df(df, results):
    # row[1] is the data row[0] is the columns
    first_entry_columns = next((row[0] for row in results if notnull(row[0])), [])
    columns = [col.name for col in first_entry_columns]
    elements = len(columns)
    data = [row[1] if notnull(row[1]) else (None,) * elements for row in results]
    dc_df = pd.DataFrame(data=data, columns=columns, index=df.index)
    dc_df = make_str_floor_to_int_number(dc_df, "id")
    return dc_df
