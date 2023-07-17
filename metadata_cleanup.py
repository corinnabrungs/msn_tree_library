import logging

import pandas as pd
from tqdm import tqdm

import lotus_client
import unichem_client
from broadinstitute_client import broad_list_search
from chembl_client import chembl_search_id_and_inchikey
from drug_utils import get_clinical_phase_description
from drugbank_client import drugbank_search_add_columns
from drugcentral_client import drugcentral_search
from lotus_client import apply_search_on_split_inchikey
from meta_constants import MetaColumns
from npatlas_client import search_np_atlas
from pandas_utils import read_dataframe, add_filename_suffix, get_parquet_file, save_dataframe, remove_empty_strings, \
    update_dataframes, remove_line_breaks
from pubchem_client import pubchem_search_structure_by_name, pubchem_search_by_structure, \
    pubchem_search_structure_by_cid, pubchem_search_parent
from rdkit_mol_identifiers import clean_structure_add_mol_id_columns, ensure_smiles_column
from structure_classifier_client import apply_np_classifier, apply_classyfire
from synonyms import ensure_synonyms_column, extract_synonym_ids

tqdm.pandas()
logging.basicConfig(format="%(asctime)s - %(message)s", level=logging.DEBUG)


def extract_prepare_input_data(metadata_file, lib_id, plate_id_header="plate_id", well_header="well_location",
                               use_cached_parquet_file: bool = True):
    """
    Load prepared metadata file from parquet or if not present prepare metadata file. Creates a row_id column to later
    align results with
    :param metadata_file: file to load
    :param lib_id: library id unique description
    :param plate_id_header: header for plate_id
    :param well_header: header for well location
    :param use_cached_parquet_file: True, try to load a parquet file with intermediate results. Otherwise overwrite this file
    :return:
    """
    out_file = get_parquet_file(metadata_file)
    if use_cached_parquet_file:
        try:
            df = read_dataframe(out_file)
            logging.info("Loaded prepared metadata file from " + out_file)
            return df
        except:
            logging.info("No cached parquet file found - will load original file")
            pass

    # load original data
    logging.info("Preparing metadata file " + metadata_file)
    df = read_dataframe(metadata_file)
    create_unique_sample_id_column(df, lib_id, plate_id_header, well_header)
    # needed as we are going to duplicate some rows if there is conflicting information, e.g., the structure in PubChem
    df["row_id"] = df.reset_index().index
    # ensure compound_name is saved to input_name if this is not defined yet
    if MetaColumns.input_name not in df.columns and MetaColumns.compound_name in df.columns:
        df[MetaColumns.input_name] = df[MetaColumns.compound_name]

    ensure_smiles_column(df)

    # apply this here to ensure that structures given as inchi are also present as smiles
    # get mol from smiles or inchi
    # calculate all identifiers from mol - monoisotopic_mass, ...
    df = clean_structure_add_mol_id_columns(df, drop_mol=True)

    df.loc[
        df[MetaColumns.smiles].notnull() | df[MetaColumns.inchi].notnull(), MetaColumns.structure_source] = "input"

    df = ensure_synonyms_column(df)
    logging.info("Writing prepared metadata file to " + out_file)
    df.to_parquet(out_file)
    return df


def save_intermediate_parquet(df, metadata_file):
    """
    Overwrite parquet file
    """
    df = df.drop(columns=["pubchem"], errors="ignore")
    save_dataframe(df, get_parquet_file(metadata_file))


def save_results(df, metadata_file):
    """
    Overwrite parquet file and create cleaned metadata file as tsv or csv (depending on input)
    """
    df = df.apply(lambda v: remove_line_breaks(v))
    save_dataframe(df, add_filename_suffix(metadata_file, "cleaned"))
    save_intermediate_parquet(df, metadata_file)


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
                df[MetaColumns.unique_sample_id] = ["{}_{}_{}_id".format(lib_id, plate, well) for plate, well in
                                                    zip(df[plate_id_header], df[well_header])]
            else:
                df[MetaColumns.unique_sample_id] = ["{}_{}_id".format(lib_id, well) for well in df[well_header]]
    except:
        logging.info(
            "No well location and plate id found to construct unique ID which is needed for library generation")
        pass


def drop_duplicates_by_structure_rowid_reset_index(df):
    try:
        id_columns = ["inchikey", "row_id"]
        df = df.drop_duplicates(id_columns, keep="first").sort_index().reset_index(drop=True)
    except:
        pass
    return df


def search_all_unichem_xrefs(df, metadata_file):
    """
    Search unichem, save to file, add id columns to dataframe
    """
    unichem_df = unichem_client.search_all_xrefs(df)
    if unichem_df is None or len(unichem_df) == 0:
        return df

    unichem_client.save_unichem_df(metadata_file, unichem_df)
    df = unichem_client.extract_ids_to_columns(unichem_df, df)
    return df


def cleanup_file(metadata_file, lib_id, plate_id_header="plate_id", well_header="well_location",
                 use_cached_parquet_file: bool = True,
                 query_pubchem_by_cid: bool = True,
                 query_pubchem_by_name: bool = True, calc_identifiers: bool = True, query_unichem: bool = True,
                 query_pubchem_by_structure: bool = True, query_chembl: bool = True, query_npclassifier: bool = True,
                 query_classyfire: bool = True, query_npatlas: bool = True, query_broad_list: bool = False,
                 query_drugbank_list: bool = False, query_drugcentral: bool = False, query_lotus: bool = False):
    logging.info("Will run on %s", metadata_file)
    df = extract_prepare_input_data(metadata_file, lib_id, plate_id_header, well_header, use_cached_parquet_file)

    if query_pubchem_by_cid:
        df = pubchem_search_parent(df, apply_structures=True)

    save_intermediate_parquet(df, metadata_file)

    # Query pubchem by name and CAS
    if query_pubchem_by_name:
        df = pubchem_search_structure_by_name(df)

    save_intermediate_parquet(df, metadata_file)

    # get mol from smiles or inchi
    # calculate all identifiers from mol - monoisotopic_mass, ...
    if calc_identifiers:
        df = clean_structure_add_mol_id_columns(df, drop_mol=True)

    if query_pubchem_by_structure:
        df = pubchem_search_by_structure(df)
        save_intermediate_parquet(df, metadata_file)

    if query_pubchem_by_cid:
        df = pubchem_search_parent(df, apply_structures=True)
        if calc_identifiers:
            df = clean_structure_add_mol_id_columns(df, drop_mol=True)

        save_intermediate_parquet(df, metadata_file)

    # drop duplicates because PubChem name search might generate new rows for conflicting smiles structures
    df = drop_duplicates_by_structure_rowid_reset_index(df)
    save_intermediate_parquet(df, metadata_file)

    # add new columns for cross references to other databases
    if query_unichem:
        search_all_unichem_xrefs(df, metadata_file)

    save_intermediate_parquet(df, metadata_file)

    # extract ids like the UNII, ...
    df = ensure_synonyms_column(df)
    df = extract_synonym_ids(df)

    if query_chembl:
        df = chembl_search_id_and_inchikey(df)
        save_intermediate_parquet(df, metadata_file)

    if query_broad_list:
        df = broad_list_search(df)
        save_intermediate_parquet(df, metadata_file)

    if query_drugbank_list:
        df = drugbank_search_add_columns(df)
        save_intermediate_parquet(df, metadata_file)

    if query_drugcentral:
        df = drugcentral_search(df)
        save_intermediate_parquet(df, metadata_file)

    # Converting numbers back to phase X, launched or preclinic
    df["clinical_phase_description"] = [get_clinical_phase_description(number) for number in df["clinical_phase"]]
    # drop mol
    df = df.drop(columns=['mol', 'pubchem'], errors="ignore")
    df["none"] = df.isnull().sum(axis=1)
    try:
        df = df.sort_values(by="none", ascending=True).drop_duplicates(
            ["inchikey", "monoisotopic_mass", "row_id"], keep="first").sort_index()
    except:
        pass

    # GNPS cached version
    if query_npclassifier:
        result = apply_np_classifier(df)
        df = update_dataframes(result, df)
        save_intermediate_parquet(df, metadata_file)

    # GNPS cached version
    if query_classyfire:
        result = apply_classyfire(df)
        df = update_dataframes(result, df)
        save_intermediate_parquet(df, metadata_file)

    if query_npatlas:
        result = search_np_atlas(df)
        df = update_dataframes(result, df)
        save_intermediate_parquet(df, metadata_file)

    if query_lotus:
        df = lotus_client.apply_search_on_split_inchikey(df)

    # export metadata file
    save_results(df, metadata_file)


if __name__ == "__main__":
    cleanup_file(r"examples\test_metadata_small.tsv", lib_id="test",
                 use_cached_parquet_file=True,
                 query_pubchem_by_name=True,
                 # need local files
                 query_broad_list=True, query_drugbank_list=True, query_drugcentral=True, query_lotus=True)
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
