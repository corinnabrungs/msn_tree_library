import logging
from dataclasses import dataclass
from datetime import timedelta
from enum import Enum, auto
from strenum import StrEnum

import pandas as pd

from date_utils import create_expired_entries_dataframe, iso_datetime_now
from meta_constants import MetaColumns
from pandas_utils import (
    get_first_value_or_else,
    update_dataframes,
    save_dataframe,
    read_dataframe,
    notnull,
)
from tqdm import tqdm

tqdm.pandas()
from rest_utils import get_json_response_with_headers
from joblib import Memory

memory = Memory("memcache")

# links to the unichem interface, insert a UCI
UNI_CHEM_COMPOUND_URL = (
    "https://www.ebi.ac.uk/unichem/compoundsources?type=uci&compound={}"
)
#
UNI_CHEM_URL = "https://www.ebi.ac.uk/unichem/api/v1/compounds"
UNI_CHEM_SUFFIX = "_unichem"


class Columns(StrEnum):
    shortName = auto()
    url = auto()
    compoundId = auto()
    inchikey = auto()


@dataclass
class Source:
    shortname: str
    column_name: str = None

    def __post_init__(self):
        if self.column_name is None:
            self.column_name = self.shortname + "_id"


class Sources(Source, Enum):
    pubchem_source = ("pubchem", "pubchem_cid")
    chembl_source = "chembl"
    drugbank_source = "drugbank"
    hmdb_source = "hmdb"
    zinc_source = "zinc"
    surechembl_source = ("surechembl", "schembl_id")
    unichem_source = "unichem"
    nmrshiftdb2_source = "nmrshiftdb2"
    kegg_ligand_source = "kegg_ligand"
    drugcentral_source = "drugcentral"
    comptox_source = "comptox"
    mcule_source = "mcule"
    chebi_source = "chebi"
    lincs_source = "lincs"
    gtopdb_source = "gtopdb"  # guide to pharmacology
    rxnorm_source = "rxnorm"
    probes_and_drugs_source = "probes_and_drugs"
    fdasis_source = ("fdasis", "unii")


def search_unichem_xref(structure: str, search_type="inchikey") -> None | dict:
    try:
        return _search_unichem_xref(structure, search_type)
    except:
        logging.info(f"Failed unichem for {structure} as {search_type}")
        return None


@memory.cache
def _search_unichem_xref(structure: str, search_type="inchikey") -> None | dict:
    if structure is None or len(structure) == 0:
        return None

    body = {"type": search_type, "compound": structure}
    headers = {"Content-Type": "application/json"}
    result = get_json_response_with_headers(
        UNI_CHEM_URL, headers, body, True, timeout=5
    )

    if result is None:
        raise Exception("unichem service failed")

    if result is not None and "compounds" in result and len(result["compounds"]) > 0:
        for comp in result["compounds"]:
            uci_str = str(comp["uci"])
            uci_source = {
                "compoundId": uci_str,
                "shortName": "unichem",
                "url": UNI_CHEM_COMPOUND_URL.format(uci_str),
                "longName": "UniChem",
            }
            comp["sources"].append(uci_source)

    return result


def search_unichem_xref_for_row(row) -> pd.DataFrame:
    """
    Maybe multiple matches as rows
    :param row:
    :return:
    """
    # try inchikey - otherwise smiles
    result = search_unichem_xref(row["inchikey"], "inchikey")

    if result is not None and "compounds" in result and len(result["compounds"]) > 0:
        dataframes = [pd.DataFrame(comp["sources"]) for comp in result["compounds"]]
        result_df = pd.concat(dataframes)
        result_df["inchikey"] = row["inchikey"]
        return result_df[["inchikey", "compoundId", "shortName", "url"]]
    else:
        # return empty dataframe
        return pd.DataFrame(
            {"inchikey": [], "compoundId": [], "shortName": [], "url": []}
        )


def search_all_xrefs(
    df: pd.DataFrame,
    metadata_file=None,
    refresh_expired_entries_after: timedelta = timedelta(days=90),
) -> pd.DataFrame:
    """

    :param df: requires columns inchikey
    :return: a new data frame with all xrefs and their inchikey (input)
    """
    logging.info("Running uni chem search for cross references")

    # only work on expired elements
    # define which rows are old or were not searched before
    filtered = create_expired_entries_dataframe(
        df, MetaColumns.date_unichem_search, refresh_expired_entries_after
    )

    if len(filtered) == 0:  # no need to update
        return pd.DataFrame()  # return empty dataframe

    unichem_results = pd.concat(
        filtered.progress_apply(
            lambda row: search_unichem_xref_for_row(row), axis=1
        ).array
    )
    # unichem_results = unichem_results.drop_duplicates().reset_index(drop=True)
    if unichem_results is None or len(unichem_results) == 0:
        return df

    if metadata_file is not None:
        save_unichem_df(metadata_file, unichem_results)

    return extract_ids_to_columns(unichem_results, filtered)


def extract_ids_to_columns(
    unichem_df: pd.DataFrame, target_df: pd.DataFrame
) -> pd.DataFrame:
    """

    :param unichem_df: contains all unichem results with inchikey and ids
    :param target_df: metadata contains inchikey and this function adds new columns
    :return: the metadata target_df
    """
    results = target_df[[MetaColumns.inchikey]].copy()
    for source in Sources:
        col_name = source.column_name
        results[col_name] = [
            get_compoundid(unichem_df, inchikey, source)
            for inchikey in results[MetaColumns.inchikey]
        ]

    unichem_url = Sources.unichem_source.shortname + "_url"
    results[unichem_url] = [
        get_unichem_url(unichem_df, inchikey)
        for inchikey in results[MetaColumns.inchikey]
    ]

    date_now = iso_datetime_now()
    target_df[MetaColumns.date_unichem_search] = [
        date_now if notnull(res) else None for res in results[MetaColumns.unichem_id]
    ]

    return update_dataframes(results, target_df)


def get_unichem_url(unichem_df: pd.DataFrame, inchikey: str) -> str | None:
    return get_source_url(unichem_df, inchikey, Sources.unichem_source)


def get_compoundid(
    unichem_df: pd.DataFrame, inchikey: str, source: Source
) -> str | None:
    return _get_first_value(unichem_df, inchikey, source, Columns.compoundId)


def get_source_url(
    unichem_df: pd.DataFrame, inchikey: str, source: Source
) -> str | None:
    return _get_first_value(unichem_df, inchikey, source, Columns.url)


def _get_first_value(
    unichem_df: pd.DataFrame, inchikey: str, source: Source, extract_column: Columns
) -> str | None:
    df = unichem_df.query(f"inchikey=='{inchikey}' & shortName=='{source.shortname}'")
    if len(df) == 0:
        return None

    return get_first_value_or_else(df, extract_column)


def save_unichem_df(base_file, df: pd.DataFrame):
    from pandas_utils import add_filename_suffix

    parquet_file = add_filename_suffix(base_file, "unichem", ".parquet.gzip")
    # try read old results and only replace new ones
    try:
        old = read_dataframe(parquet_file).set_index(MetaColumns.unichem_id)
        df = df.set_index(MetaColumns.unichem_id)
        df = update_dataframes(df, old).reset_index()
    except:
        logging.info("No unichem file found - creating a new one at " + parquet_file)

    save_dataframe(df, add_filename_suffix(base_file, "unichem", ".csv"))
    save_dataframe(df, parquet_file)
