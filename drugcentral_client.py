import logging
import pandas as pd
from tqdm import tqdm
from pandas_utils import notnull

import drugcentral_postgresql_query as drugcentral_query
from rdkit_mol_identifiers import split_inchikey

tqdm.pandas()


def drugcentral_search(df):
    logging.info("Search drugcentral by external identifier or inchikey")
    prefix = "drugcentral_"
    if "split_inchikey" not in df and "inchikey" in df:
        df["split_inchikey"] = [split_inchikey(inchikey) for inchikey in df['inchikey']]
    try:
        drugcentral_query.connect()
        logging.info("Searching in DrugCentral")
        results = df.progress_apply(lambda row: drugcentral_query.drugcentral_for_row(row), axis=1)
        # results = [drugcentral_query.drugcentral_for_row(row) for _,row in df.iterrows()]
        # results = [drugcentral_query.drugcentral_postgresql(inchikey, split_inchikey) for inchikey, split_inchikey in
        #          tqdm(zip(df["inchikey"], df["split_inchikey"]))]
        logging.info("DrugCentral search done")

        # row[1] is the data row[0] is the columns
        first_entry_columns = next((row[0] for row in results if notnull(row[0])), [])
        columns = [col.name for col in first_entry_columns]
        elements = len(columns)
        data = [row[1] if notnull(row[1]) else (None,) * elements for row in results]
        dc_df = pd.DataFrame(data=data, columns=columns, index=df.index)
        dc_df = dc_df.add_prefix(prefix)
        df = pd.concat([df, dc_df], axis=1)

        if "drugcentral_administration" in df.columns:
            df["drugcentral_administration_number"] = [4 if notnull(status) else None for status in
                                                       df["drugcentral_administration"]]
        else:
            df["drugcentral_administration_number"] = None
        df["clinical_phase"] = df[['clinical_phase', 'drugcentral_administration_number']].max(axis=1)
        if "drugbank_approved" in df.columns:
            df["any_phase"] = df["drugbank_approved"].notna() | (df["clinical_phase"] > 0)
        else:
            df["any_phase"] = df["clinical_phase"] > 0

        return df
    finally:
        drugcentral_query.deconnect()
