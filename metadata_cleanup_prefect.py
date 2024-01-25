import logging

from tqdm import tqdm
import argparse

import lotus_client
from broadinstitute_client import broad_list_search
from chembl_client import chembl_search_id_and_inchikey
from drug_utils import get_clinical_phase_description, map_clinical_phase_to_number
from drugbank_client import drugbank_search_add_columns
from drugcentral_client import drugcentral_search
from metadata_cleanup import (
    extract_prepare_input_data,
    save_results,
    drop_duplicates_by_structure_rowid_reset_index,
    search_all_unichem_xrefs,
    save_intermediate_parquet,
)
from npatlas_client import search_np_atlas
from pandas_utils import update_dataframes, create_missing_columns
from pubchem_client import (
    pubchem_search_structure_by_name,
    pubchem_search_by_structure,
    pubchem_search_structure_by_cid,
    pubchem_search_parent,
)
from rdkit_mol_identifiers import clean_structure_add_mol_id_columns
from synonyms import ensure_synonyms_column, extract_synonym_ids
from structure_classifier_client import apply_np_classifier, apply_classyfire
from prefect import flow, task

tqdm.pandas()
logging.basicConfig(format="%(asctime)s - %(message)s", level=logging.DEBUG)


@task(name="Extract prepare data")
def extract_prepare_input_data_prefect(
    metadata_file,
    lib_id,
    plate_id_header="plate_id",
    well_header="well_location",
    use_cached_parquet_file: bool = True,
):
    return extract_prepare_input_data(
        metadata_file, lib_id, plate_id_header, well_header, use_cached_parquet_file
    )


@task(name="save results")
def save_results_prefect(df, metadata_file):
    save_results(df, metadata_file)


@task(name="save intermediate results")
def save_intermediate_parquet_prefect(df, metadata_file):
    save_intermediate_parquet(df, metadata_file)


@task(name="PubChem parent by CID")
def pubchem_search_parent_by_cid_prefect(df, apply_structures: bool):
    return pubchem_search_parent(df, apply_structures)


@task(name="PubChem by name")
def pubchem_search_structure_by_name_prefect(df):
    return pubchem_search_structure_by_name(df)


@task(name="Clean structures")
def clean_structure_add_mol_id_columns_prefect(df):
    return clean_structure_add_mol_id_columns(df, drop_mol=True)


@task(name="search_all_unichem_xrefs")
def search_all_unichem_xrefs_prefect(df, metadata_file):
    return search_all_unichem_xrefs(df, metadata_file)


@task(name="pubchem_search_by_structure")
def pubchem_search_by_structure_prefect(df):
    return pubchem_search_by_structure(df)


@task(name="chembl_search_id_and_inchikey")
def chembl_search_id_and_inchikey_prefect(df):
    return chembl_search_id_and_inchikey(df)


@task(name="broad_list_search")
def broad_list_search_prefect(df):
    return broad_list_search(df)


@task(name="drugbank_search_add_columns")
def drugbank_search_add_columns_prefect(df):
    return drugbank_search_add_columns(df)


@task(name="drugcentral_search")
def drugcentral_search_prefect(df):
    return drugcentral_search(df)


@task(name="apply_np_classifier")
def apply_np_classifier_prefect(df):
    return apply_np_classifier(df)


@task(name="apply_classyfire")
def apply_classyfire_prefect(df):
    return apply_classyfire(df)


@task(name="apply_npatlas")
def search_np_atlas_prefect(df):
    return search_np_atlas(df)


@task(name="search lotus compound taxon links")
def search_lotus_prefect(df):
    return lotus_client.apply_search_on_split_inchikey(df)


@flow(name="Add lotus", version="0.1.0", flow_run_name="{lib_id}:{metadata_file}")
def add_lotus_flow(
    metadata_file, lib_id, plate_id_header="plate_id", well_header="well_location"
):
    logging.info("Will run on %s", metadata_file)
    df = extract_prepare_input_data_prefect(
        metadata_file, lib_id, plate_id_header, well_header
    )
    df = search_lotus_prefect(df)
    save_results_prefect(df, metadata_file)


@flow(
    name="Metadata cleanup", version="0.1.0", flow_run_name="{lib_id}:{metadata_file}"
)
def cleanup_file(
    metadata_file,
    lib_id,
    plate_id_header="plate_id",
    well_header="well_location",
    use_cached_parquet_file: bool = True,
    query_pubchem_by_cid: bool = True,
    query_pubchem_by_name: bool = True,
    calc_identifiers: bool = True,
    query_unichem: bool = True,
    query_pubchem_by_structure: bool = True,
    query_chembl: bool = True,
    query_npclassifier: bool = True,
    query_classyfire: bool = True,
    query_npatlas: bool = True,
    query_broad_list: bool = False,
    query_drugbank_list: bool = False,
    query_drugcentral: bool = False,
    query_lotus: bool = False,
):
    logging.info("Will run on %s", metadata_file)
    df = extract_prepare_input_data_prefect(
        metadata_file, lib_id, plate_id_header, well_header, use_cached_parquet_file
    )

    if query_pubchem_by_cid:
        df = pubchem_search_parent_by_cid_prefect(df, apply_structures=True)

    # Query pubchem by name and CAS
    if query_pubchem_by_name:
        df = pubchem_search_structure_by_name_prefect(df)

    save_intermediate_parquet_prefect(df, metadata_file)

    # get mol from smiles or inchi
    # calculate all identifiers from mol - monoisotopic_mass, ...
    if calc_identifiers:
        df = clean_structure_add_mol_id_columns_prefect(df)

    if query_pubchem_by_structure:
        df = pubchem_search_by_structure_prefect(df)
        save_intermediate_parquet_prefect(df, metadata_file)

    if query_pubchem_by_cid:
        df = pubchem_search_parent_by_cid_prefect(df, apply_structures=True)
        if calc_identifiers:
            df = clean_structure_add_mol_id_columns_prefect(df)
        save_intermediate_parquet_prefect(df, metadata_file)

    # drop duplicates because PubChem name search might generate new rows for conflicting smiles structures
    df = drop_duplicates_by_structure_rowid_reset_index(df)
    save_intermediate_parquet_prefect(df, metadata_file)

    # structures are now fetched. Run things in parallel
    # run in parallel
    tasks = []

    # GNPS cached version
    if query_npclassifier:
        tasks.append(apply_np_classifier_prefect.submit(df))

    # GNPS cached version
    if query_classyfire:
        tasks.append(apply_classyfire_prefect.submit(df))

    if query_npatlas:
        tasks.append(search_np_atlas_prefect.submit(df))

    if query_lotus:
        tasks.append(search_lotus_prefect.submit(df))

    # add new columns for cross references to other databases
    if query_unichem:
        # xrefs are needed for other steps so run sequential here
        df = search_all_unichem_xrefs_prefect(df, metadata_file)
        save_intermediate_parquet_prefect(df, metadata_file)

    # extract ids like the UNII, ...
    df = ensure_synonyms_column(df)
    df = extract_synonym_ids(df)

    if query_chembl:
        df = chembl_search_id_and_inchikey_prefect(df)
        save_intermediate_parquet_prefect(df, metadata_file)

    if query_broad_list:
        df = broad_list_search_prefect(df)

    if query_drugbank_list:
        df = drugbank_search_add_columns_prefect(df)

    if query_drugcentral:
        df = drugcentral_search_prefect(df)

    # Converting numbers
    df = create_missing_columns(df, ["clinical_phase", "Clinical Information"])
    df["clinical_phase"] = [
        map_clinical_phase_to_number(phase) for phase in df["clinical_phase"]
    ]
    df["Clinical Information_clinical_phase"] = [
        map_clinical_phase_to_number(phase) for phase in df["Clinical Information"]
    ]

    # Getting highest number
    df["clinical_phase"] = df[
        df.columns[df.columns.str.endswith("clinical_phase")]
    ].max(axis=1)

    # Converting numbers back to phase X, launched or preclinic
    df["clinical_phase_description"] = [
        get_clinical_phase_description(number) for number in df["clinical_phase"]
    ]

    df["any_phase"] = df["clinical_phase"] > 0
    # drop mol
    df = df.drop(columns=["mol", "pubchem"], errors="ignore")
    df["none"] = df.isnull().sum(axis=1)
    try:
        df = (
            df.sort_values(by="none", ascending=True)
            .drop_duplicates(["inchikey", "monoisotopic_mass", "row_id"], keep="first")
            .sort_index()
        )
    except:
        pass

    # # wait for all tasks in parallel to finish
    for task in tasks:
        # get results and merge into
        result_df = task.result()
        df = update_dataframes(result_df, df)

    # result_dfs = [task.result() for task in tasks]
    # df = df.join(result_dfs, sort=False)

    # export metadata file
    save_results_prefect(df.copy(), metadata_file)


def full_cleanup_file(metadata_file, lib_id, use_cached_parquet_file: bool = True):
    try:
        cleanup_file(
            metadata_file,
            lib_id,
            use_cached_parquet_file=use_cached_parquet_file,
            query_pubchem_by_name=True,
            # need local files
            query_broad_list=True,
            query_drugbank_list=False,
            query_drugcentral=True,
            query_lotus=True,
        )
    except:
        logging.exception("Exception in flow")
        # exit(1)




if __name__ == "__main__":
    # full_cleanup_file(r"examples\test_metadata.tsv", lib_id="test")
    full_cleanup_file(r"examples\test_metadata_small.tsv", lib_id="test")
    exit(0)

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "input_file",
        type=str,
        help="the input metadata file, with some of these columns inchikey, smiles, inchi, compound_name, synonyms, ",
    )
    parser.add_argument(
        "-l", "--lib_id", type=str, help="library id is added to some columns"
    )
    args = parser.parse_args()

    full_cleanup_file(args.input_file, args.lib_id)
    exit(0)


    # cleanup_file("data\mce_library.tsv", id_columns=['Product Name', 'lib_plate_well', "inchikey"], query_pubchem=True, query_broad_list=True, query_drugbank_list=True,
    #                  query_drugcentral=True)
    # cleanup_file("data\gnpslib\gnps_library_small.csv", id_columns=['gnps_libid', "inchikey"], query_pubchem=True,
    #              query_broad_list=True, query_drugbank_list=True, query_drugcentral=True)
    # cleanup_file("data\gnpslib\gnps_library.csv", id_columns=['gnps_libid', "inchikey"], query_pubchem=True,
    #              query_broad_list=True, query_drugbank_list=True, query_drugcentral=True)
    # cleanup_file("data\mce_library_add_compounds.tsv", id_columns=['Product Name', 'lib_plate_well', "inchikey"], query_pubchem=True, query_broad_list=True, query_drugbank_list=True,
    #                  query_drugcentral=True)
    # cleanup_file(r"data/nih/nih_library_test.tsv", id_columns=['Product Name', 'lib_plate_well', "smiles"], query_pubchem=True,  pubchem_search=True, query_broad_list=True, query_drugbank_list=True,
    #                  query_drugcentral=True)
