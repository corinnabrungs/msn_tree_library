import logging
from dataclasses import dataclass
from enum import Enum, auto
from strenum import StrEnum

import pandas as pd

from rest_utils import get_json_response_with_headers

# links to the unichem interface, insert a UCI
UNI_CHEM_COMPOUND_URL = "https://www.ebi.ac.uk/unichem/compoundsources?type=uci&compound={}"
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
    surechembl_source = ("surechembl", "schembl")
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
    if structure is None or len(structure) == 0:
        return None

    body = {
        "type": search_type,
        "compound": structure
    }
    headers = {"Content-Type": "application/json"}
    result = get_json_response_with_headers(UNI_CHEM_URL, headers, body, True)

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
            {
                "inchikey": [],
                "compoundId": [],
                "shortName": [],
                "url": []
            }
        )


def search_all_xrefs(df) -> pd.DataFrame:
    """

    :param df: requires columns inchikey
    :return: a new data frame with all xrefs and their inchikey (input)
    """
    logging.info("Running uni chem search for cross references")
    unichem_results = pd.concat(df.apply(lambda row: search_unichem_xref_for_row(row), axis=1).array)
    return unichem_results


def extract_ids_to_columns(unichem_df: pd.DataFrame, target_df: pd.DataFrame) -> pd.DataFrame:
    """

    :param unichem_df: contains all unichem results with inchikey and ids
    :param target_df: metadata contains inchikey and this function adds new columns
    :return: the metadata target_df
    """
    for source in Sources:
        col_name = source.column_name
        target_df[col_name] = [get_compoundid(unichem_df, inchikey, source) for inchikey in target_df["inchikey"]]

    unichem_url = Sources.unichem_source.shortname + "_url"
    target_df[unichem_url] = [get_unichem_url(unichem_df, inchikey) for inchikey in target_df["inchikey"]]

    return target_df


def get_unichem_url(unichem_df: pd.DataFrame, inchikey: str) -> str | None:
    return get_source_url(unichem_df, inchikey, Sources.unichem_source)


def get_compoundid(unichem_df: pd.DataFrame, inchikey: str, source: Source) -> str | None:
    return _get_first_value(unichem_df, inchikey, source, Columns.compoundId)


def get_source_url(unichem_df: pd.DataFrame, inchikey: str, source: Source) -> str | None:
    return _get_first_value(unichem_df, inchikey, source, Columns.url)


def _get_first_value(unichem_df: pd.DataFrame, inchikey: str, source: Source, extract_column: Columns) -> str | None:
    df = unichem_df.query(f"inchikey=='{inchikey}' & shortName=='{source.shortname}'")
    if len(df) == 0:
        return None

    return df.at[df.index[0], extract_column]


def save_unichem_df(base_file, df: pd.DataFrame):
    from pandas_utils import add_filename_suffix
    df.to_csv(add_filename_suffix(base_file, "unichem", ".csv"), index=False)
    df.to_parquet(add_filename_suffix(base_file, "unichem", ".parquet.gzip"), compression='gzip')
