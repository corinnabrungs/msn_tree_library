import pandas as pd

from date_utils import iso_datetime_now
from meta_constants import MetaColumns
from pandas_utils import read_dataframe, update_dataframes

lotus_prefix = "lotus_"


def read_lotus_dataframe(use_unique_split_inchikeys: bool = True) -> pd.DataFrame:
    if use_unique_split_inchikeys:
        return read_dataframe("data/lotus_unique_split_inchikey_download.parquet")
    else:
        return read_dataframe("data/lotus_download.parquet")


def apply_search_on_split_inchikey(df: pd.DataFrame) -> pd.DataFrame:
    lotus = read_lotus_dataframe(use_unique_split_inchikeys=True).drop(columns=["inchikey"])
    lotus = lotus.add_prefix(lotus_prefix)
    lotus = lotus.rename(columns={lotus_prefix + MetaColumns.split_inchikey: MetaColumns.split_inchikey})

    # create df with only split_inchikey and merge in
    results = df[[MetaColumns.split_inchikey]].copy()
    results = results.merge(lotus, how="left", sort=False, on=MetaColumns.split_inchikey)

    # refresh date
    results.loc[results[lotus_prefix + "taxon"].notnull(), MetaColumns.date_wikidata_lotus_search] = iso_datetime_now()
    return update_dataframes(results, df)
