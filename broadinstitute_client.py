import pandas as pd
import logging
from tqdm import tqdm

import pandas_utils
from date_utils import iso_datetime_now
from meta_constants import MetaColumns
from pandas_utils import left_merge_retain_index, add_column_prefix, update_dataframes, create_missing_columns

tqdm.pandas()


def broad_list_search(df):
    from rdkit_mol_identifiers import split_inchikey
    from drug_utils import map_clinical_phase_to_number

    if MetaColumns.split_inchikey not in df and MetaColumns.inchikey in df:
        df[MetaColumns.split_inchikey] = [split_inchikey(inchikey) for inchikey in df[MetaColumns.inchikey]]
    # only merge on id column where id is notnull
    results = df[[MetaColumns.split_inchikey]][df[MetaColumns.split_inchikey].notnull()].copy()
    if len(results) == 0:
        return df

    logging.info("Search broad institute list of drugs by first block of inchikey")
    # download from: https://clue.io/repurposing#download-data
    prefix = "broad_"
    broad_df = pd.read_csv("data/broad_institute_drug_list.csv")

    broad_df[MetaColumns.split_inchikey] = [split_inchikey(inchikey) for inchikey in broad_df["InChIKey"]]
    broad_df = broad_df.drop(columns=["InChIKey"])

    results = left_merge_retain_index(results, broad_df, on=MetaColumns.split_inchikey)
    results = add_column_prefix(results, prefix, columns_to_keep=MetaColumns.split_inchikey)

    # converting the clinical phases (from broad institute, chembl, provider, or else) to numbers (remove phase,
    # preclinic (as 0.5), or launched)
    # TODO do this somewhere esle
    results = create_missing_columns(results, ['broad_clinical_phase', 'clinical_phase', 'Clinical Information'])

    results["broad_clinical_phase"] = [map_clinical_phase_to_number(phase) for phase in
                                       results["broad_clinical_phase"]]
    results["clinical_phase"] = [map_clinical_phase_to_number(phase) for phase in results["clinical_phase"]]
    results["Clinical Information"] = [map_clinical_phase_to_number(phase) for phase in results["Clinical Information"]]

    # Comparing the clinical phases, only store the highest number in clinical phase column
    results["clinical_phase"] = results[['broad_clinical_phase', 'clinical_phase', 'Clinical Information']].max(axis=1)

    results = results.drop(columns=["broad_clinical_phase", "Clinical Information"], errors="ignore")
    results[MetaColumns.date_broad_drug_list] = iso_datetime_now()
    return update_dataframes(results, df)
