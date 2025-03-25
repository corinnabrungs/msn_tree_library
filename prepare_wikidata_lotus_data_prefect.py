import lotus_client
import pandas_utils
import rdkit_mol_identifiers
from meta_constants import MetaColumns
from metadata_cleanup_prefect import clean_structure_add_mol_id_columns_prefect
from pandas_utils import (
    notnull,
    isnull,
    divide_chunks,
    groupby_join_unique_values,
    isnull_or_empty,
)
import pandas as pd
import numpy as np
from SPARQLWrapper import SPARQLWrapper, JSON
import sys
from rdkit_mol_identifiers import split_inchikey, clean_structure_add_mol_id_columns

wikidata_sparql_url = "https://query.wikidata.org/sparql"
from prefect import task, flow


def ensure_identifier(item, identifier: str = "wd"):
    """

    :param item: may be url or item or already identified: http://www.wikidata.org/entity/Q193572 or Q193572
    :param identifier:
    :return: identifier:item e.g. wd:Q193572
    """
    if isnull_or_empty(item):
        return None
    item = str(item)

    prefix = identifier + ":"
    if item.startswith("http"):
        o = item.split("/")
        if len(o) > 1:
            return prefix + o[-1]

    if item.startswith(prefix):
        return item

    return prefix + item


def _create_input_list(taxon_list, identifier: str = "wd") -> str:
    items = [
        ensure_identifier(item, identifier)
        for item in taxon_list
        if notnull(item) and len(str(item)) > 0
    ]
    return " ".join(items)


def _get_lotus_compound_taxon_relations_sparql():
    return """
    #title: Which are the available referenced structure-organism pairs on Wikidata?
    SELECT DISTINCT ?structure ?inchikey ?taxon ?reference WHERE {
      ?structure wdt:P235 ?inchikey;                 # get the inchikey
        p:P703[                                      # statement found in taxon
         ps:P703 ?taxon;                             # get the taxon
         (prov:wasDerivedFrom/pr:P248) ?reference ]. # get the reference
    }
    """


def _get_reference_info_sparql(references):
    input_list = _create_input_list(references, "wd")
    return """
    #title: get NCBI ids from input taxon and parent
    SELECT DISTINCT ?reference ?doi
    WHERE {
      VALUES ?reference {
        # wd:Q30046  # example
        PLACEHOLDER
      }
      OPTIONAL { ?reference wdt:P356 ?doi. }
      SERVICE wikibase:label { bd:serviceParam wikibase:language "en". }
    }
    """.replace(
        "PLACEHOLDER", input_list
    )


def _get_structure_info_sparql(structures):
    input_list = _create_input_list(structures, "wd")
    return """
    #title: get NCBI ids from input taxon and parent
    SELECT DISTINCT ?structure ?structureLabel ?inchikey ?inchi ?isomeric_smiles ?canonical_smiles ?formula
    WHERE {
      VALUES ?structure {
        # wd:Q30046  # example
        PLACEHOLDER
      }
      OPTIONAL { ?structure wdt:P235 ?inchikey. }
      OPTIONAL { ?structure wdt:P2017 ?isomeric_smiles .}
      OPTIONAL { ?structure wdt:P233 ?canonical_smiles .}
      OPTIONAL { ?structure wdt:P274 ?formula .}
      OPTIONAL { ?structure wdt:P234 ?inchi .}
      SERVICE wikibase:label { bd:serviceParam wikibase:language "en". }
    }
    """.replace(
        "PLACEHOLDER", input_list
    )


def _get_parent_taxon_ncbi_sparql(taxon_list):
    input_list = _create_input_list(taxon_list, "wd")
    return """
    #title: get NCBI ids from input taxon and parent
    SELECT DISTINCT ?taxon ?taxon_name ?taxon_rankLabel ?ncbi_id ?parent_taxon ?parent_taxon_name ?parent_taxon_rankLabel ?parent_ncbi_id 
    WHERE {
      VALUES ?taxon {
        # wd:Q30046  # example
        PLACEHOLDER
      }
      OPTIONAL { ?taxon wdt:P225 ?taxon_name. }
      OPTIONAL { ?taxon wdt:P685 ?ncbi_id. }
      OPTIONAL { ?taxon wdt:P685 ?taxon_rank. }
      ?taxon wdt:P171 ?parent_taxon.
      ?parent_taxon wdt:P225 ?parent_taxon_name.
      OPTIONAL { ?parent_taxon wdt:P105 ?parent_taxon_rank. }
      OPTIONAL { ?parent_taxon wdt:P685 ?parent_ncbi_id. }
      SERVICE wikibase:label { bd:serviceParam wikibase:language "en". }
    }
    """.replace(
        "PLACEHOLDER", input_list
    )


def get_sparql_json_results(
    sparql_query: str, endpoint_url: str = "https://query.wikidata.org/sparql"
):
    user_agent = "WDQS-example Python/%s.%s" % (
        sys.version_info[0],
        sys.version_info[1],
    )
    # TODO adjust user agent; see https://foundation.wikimedia.org/wiki/Policy:Wikimedia_Foundation_User-Agent_Policy
    sparql = SPARQLWrapper(endpoint_url, agent=user_agent)
    sparql.setQuery(sparql_query)
    sparql.setReturnFormat(JSON)
    return sparql.queryAndConvert()


def extract_values(result):
    return {key: result[key]["value"] for key in result}


def load_as_dataframe(
    sparql_query: str, endpoint_url: str = "https://query.wikidata.org/sparql"
) -> pd.DataFrame:
    results = get_sparql_json_results(sparql_query, endpoint_url)["results"]["bindings"]
    new_results = [extract_values(result) for result in results]
    df = pd.DataFrame(new_results)
    label_cols = [col for col in df.columns if col.endswith("Label")]
    df = df.rename(columns={old: old[:-5] + "_label" for old in label_cols})
    return df


@task(name="Extract lotus from wikidata")
def extract_lotus_compound_taxon_relations():
    query = _get_lotus_compound_taxon_relations_sparql()
    df = load_as_dataframe(query, wikidata_sparql_url)
    df["split_inchikey"] = [split_inchikey(inchikey) for inchikey in df["inchikey"]]
    df = df.sort_values(["inchikey"])
    return df


@task(name="parent taxon")
def query_dataframe_parent_taxon_ncbi(taxon_list):
    chunks = divide_chunks(drop_duplicates(taxon_list), 250)
    queries = [_get_parent_taxon_ncbi_sparql(chunk) for chunk in chunks]
    dfs = [load_as_dataframe(query) for query in queries]
    return pd.concat(dfs, sort=False)


@task(name="structure info")
def query_dataframe_structure_info(structures):
    chunks = divide_chunks(drop_duplicates(structures), 250)
    queries = [_get_structure_info_sparql(chunk) for chunk in chunks]
    dfs = [load_as_dataframe(query) for query in queries]

    df = pd.concat(dfs, sort=False).sort_values(["inchikey"])
    df["split_inchikey"] = [split_inchikey(inchikey) for inchikey in df["inchikey"]]
    return df


@task(name="reference info")
def query_dataframe_reference_info(references):
    chunks = divide_chunks(drop_duplicates(references), 250)
    queries = [_get_reference_info_sparql(chunk) for chunk in chunks]
    dfs = [load_as_dataframe(query) for query in queries]
    return pd.concat(dfs, sort=False)


def drop_duplicates(input_list):
    return list(set(input_list))


def query_dataframe_lotus_compound_taxon_relations():
    query = _get_lotus_compound_taxon_relations_sparql()
    df = load_as_dataframe(query, wikidata_sparql_url)
    # add more information to the taxon, structure, and reference
    df = (
        df.merge(query_dataframe_parent_taxon_ncbi(df["taxon"]), on="taxon", how="left")
        .sort_values(["parent_ncbi_id"])
        .drop_duplicates(["taxon", "structure", "reference"])
    )

    df = (
        df.merge(
            query_dataframe_reference_info(df["reference"]), on="reference", how="left"
        )
        .sort_values(["doi"])
        .drop_duplicates(["taxon", "structure", "reference"])
    )

    structures_df = query_dataframe_structure_info(df["structure"]).drop(
        columns=["inchikey"]
    )
    df = (
        df.merge(structures_df, on="structure", how="left")
        .sort_values(["inchikey"])
        .drop_duplicates(["taxon", "structure", "reference"])
    )

    return df


@task(name="save lotus dump")
def save_lotus_dump(df):
    pandas_utils.save_dataframe(df, "data/lotus_download.parquet")
    pandas_utils.save_dataframe(df, "data/lotus_download.csv")


@task(name="Create unique inchikey entries")
def save_unique_inchikey_dump(df):
    df = groupby_join_unique_values(df, columns=["split_inchikey"], as_lists=False)
    pandas_utils.save_dataframe(df, "data/lotus_unique_split_inchikey_download.parquet")
    pandas_utils.save_dataframe(df, "data/lotus_unique_split_inchikey_download.csv")


@flow(name="download lotus")
def download_lotus_prefect():
    df = extract_lotus_compound_taxon_relations()
    # add more information to the taxon, structure, and reference

    parents = query_dataframe_parent_taxon_ncbi.submit(df["taxon"])
    references = query_dataframe_reference_info.submit(df["reference"])
    structures_df = query_dataframe_structure_info.submit(df["structure"])

    df = (
        df.merge(parents.result(), on="taxon", how="left")
        .sort_values(["parent_ncbi_id"])
        .drop_duplicates(["taxon", "structure", "reference"])
    )

    df = (
        df.merge(references.result(), on="reference", how="left")
        .sort_values(["doi"])
        .drop_duplicates(["taxon", "structure", "reference"])
    )

    # inchikey is already present from df
    structures_df = structures_df.result().drop(columns=["inchikey", "split_inchikey"])
    df = (
        df.merge(structures_df, on="structure", how="left")
        .sort_values(["inchikey"])
        .drop_duplicates(["taxon", "structure", "reference"])
    )

    # clean structures and drop not needed columns
    excluded_cols = [
        MetaColumns.smiles,
        MetaColumns.canonical_smiles,
        MetaColumns.isomeric_smiles,
        MetaColumns.inchi,
    ]
    columns = [col for col in df.columns if col not in excluded_cols]
    df = clean_structure_add_mol_id_columns_prefect(df)[columns]

    save_lotus_dump(df)
    save_unique_inchikey_dump(df)
    return df


if __name__ == "__main__":
    download_lotus_prefect()
