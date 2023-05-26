import logging

from tqdm import tqdm
import argparse
import unichem_client
from broadinstitute_client import broad_list_search
from chembl_client import chembl_search_id_and_inchikey
from drug_utils import get_clinical_phase_description
from drugbank_client import drugbank_search_add_columns
from drugcentral_client import drugcentral_search
from meta_constants import MetaColumns
from pandas_utils import read_dataframe, add_filename_suffix
from pubchem_client import pubchem_search_structure_by_name, pubchem_search_by_structure
from rdkit_mol_identifiers import clean_structure_add_mol_id_columns
from synonyms import ensure_synonyms_column, extract_synonym_ids
from structure_classifier_client import apply_np_classifier, apply_classyfire
from prefect import flow, task

tqdm.pandas()
logging.basicConfig(format="%(asctime)s - %(message)s", level=logging.DEBUG)

# CONSTANTS
unique_sample_id_header = "unique_sample_id"


def create_unique_sample_id_column(df, lib_id, plate_id_header, well_header):
    """
    generates a column with unique sample IDs if column with well locations and plate IDs (optional) is available.
    :param df: metadata
    :param lib_id: library id that defines the compound library
    :param plate_id_header: number or name of the plate
    :param well_header: well location of the injection
    :return: lib_id_plate_well_id
    """
    try:
        if well_header in df.columns:
            if plate_id_header in df.columns:
                df[unique_sample_id_header] = ["{}_{}_{}_id".format(lib_id, plate, well) for plate, well in
                                               zip(df[plate_id_header], df[well_header])]
            else:
                df[unique_sample_id_header] = ["{}_{}_id".format(lib_id, well) for well in df[well_header]]
    except:
        logging.info(
            "No well location and plate id found to construct unique ID which is needed for library generation")
        pass


def drop_duplicates_by_structure_rowid(df):
    try:
        id_columns = ["inchikey", "row_id"]
        df = df.drop_duplicates(id_columns, keep="first").sort_index()
    except:
        pass
    return df


@task(name="PubChem by name")
def pubchem_search_structure_by_name_prefect(df):
    return pubchem_search_structure_by_name(df)


@task(name="Clean structures")
def clean_structure_add_mol_id_columns_prefect(df):
    return clean_structure_add_mol_id_columns(df)


@task(name="search_all_unichem_xrefs")
def search_all_unichem_xrefs_prefect(df, metadata_file):
    unichem_df = unichem_client.search_all_xrefs(df)
    unichem_client.save_unichem_df(metadata_file, unichem_df)
    df = unichem_client.extract_ids_to_columns(unichem_df, df)
    return df


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


@task(name="save results")
def save_results_prefect(df, out_file):
    logging.info("Exporting to file %s", out_file)
    if out_file.endswith(".tsv"):
        df.to_csv(out_file, sep="\t", index=False)
    else:
        df.to_csv(out_file, sep=",", index=False)


@flow(name="Metadata cleanup",
      version="0.1.0")
def cleanup_file(metadata_file, lib_id, plate_id_header="plate_id", well_header="well_location",
                 query_pubchem_by_name: bool = True, calc_identifiers: bool = True, query_unichem: bool = True,
                 query_pubchem_by_structure: bool = True, query_chembl: bool = True, query_npclassifier: bool = True,
                 query_classyfire: bool = True, query_npatlas: bool = True, query_broad_list: bool = False,
                 query_drugbank_list: bool = False, query_drugcentral: bool = False):
    logging.info("Will run on %s", metadata_file)
    out_file = add_filename_suffix(metadata_file, "cleaned")

    df = read_dataframe(metadata_file)
    create_unique_sample_id_column(df, lib_id, plate_id_header, well_header)

    # needed as we are going to duplicate some rows if there is conflicting information, e.g., the structure in PubChem
    df["row_id"] = df.reset_index().index

    # ensure compound_name is saved to input_name if this is not defined yet
    if MetaColumns.input_name not in df.columns and MetaColumns.compound_name in df.columns:
        df[MetaColumns.input_name] = df[MetaColumns.compound_name]

    df = ensure_synonyms_column(df)

    # Query pubchem by name and CAS
    if query_pubchem_by_name:
        df = pubchem_search_structure_by_name_prefect(df)

    # get mol from smiles or inchi
    # calculate all identifiers from mol - exact_mass, ...
    if calc_identifiers:
        df = clean_structure_add_mol_id_columns_prefect(df)

    # drop duplicates because PubChem name search might generate new rows for conflicting smiles structures
    df = drop_duplicates_by_structure_rowid(df)

    # structures are now fetched. Run things in parallel
    # run in parallel
    tasks = []

    # GNPS cached version
    if query_npclassifier:
        tasks.append(apply_np_classifier_prefect.submit(df))

    # GNPS cached version
    if query_classyfire:
        tasks.append(apply_classyfire_prefect.submit(df))

    # add new columns for cross references to other databases
    if query_unichem:
        # xrefs are needed for other steps so run sequential here
        df = search_all_unichem_xrefs_prefect(df, metadata_file)

    if query_pubchem_by_structure:
        df = pubchem_search_by_structure_prefect(df)

    # extract ids like the UNII, ...
    df = ensure_synonyms_column(df)
    df = extract_synonym_ids(df)

    if query_chembl:
        df = chembl_search_id_and_inchikey_prefect(df)

    if query_broad_list:
        df = broad_list_search_prefect(df)

    if query_drugbank_list:
        df = drugbank_search_add_columns_prefect(df)

    if query_drugcentral:
        df = drugcentral_search_prefect(df)

    # Converting numbers back to phase X, launched or preclinic
    df["clinical_phase_description"] = [get_clinical_phase_description(number) for number in
                                        df["clinical_phase"]]
    # drop mol
    df = df.drop(columns=['mol', 'pubchem'])
    df["none"] = df.isnull().sum(axis=1)
    try:
        df = df.sort_values(by="none", ascending=True).drop_duplicates(
            ["inchikey", "exact_mass", "row_id"], keep="first").sort_index()
    except:
        pass

    # wait for all tasks in parallel to finish
    result_dfs = [task.result() for task in tasks]

    # export metadata file
    save_results_prefect(df, out_file)


def full_cleanup_file(metadata_file, lib_id, flow_run_name=None):
    try:
        cleanup_file(metadata_file, lib_id, query_pubchem_by_name=True,
                     # need local files
                     query_broad_list=True, query_drugbank_list=True, query_drugcentral=True
                     )
    except:
        logging.exception("Exception in flow")
        # exit(1)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("input_file", type=str,
                        help="the input metadata file, with some of these columns inchikey, smiles, inchi, compound_name, synonyms, ")
    parser.add_argument("-l", "--lib_id", type=str, help="library id is added to some columns")
    args = parser.parse_args()

    full_cleanup_file(args.input_file, args.lib_id)
    exit(0)

    # full_cleanup_file(r"examples\test_metadata.tsv", lib_id="test")
    # full_cleanup_file(r"examples\test_metadata_small.tsv", lib_id="test")

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
