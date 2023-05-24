import logging

from tqdm import tqdm

import unichem_client
from broadinstitute_client import broad_list_search
from chembl_client import chembl_search_id_and_inchikey
from drugbank_client import drugbank_search_add_columns
from drugcentral_client import drugcentral_search
from pandas_utils import read_dataframe, add_filename_suffix
from pubchem_client import pubchem_search_structure_by_name, pubchem_search_by_structure
from rdkit_mol_identifiers import clean_structure_add_mol_id_columns
from synonyms import ensure_synonyms_column, extract_synonym_ids
from structure_classifier_client import apply_np_classifier, apply_classyfire

tqdm.pandas()
logging.basicConfig(format="%(asctime)s - %(message)s", level=logging.DEBUG)

# CONSTANTS
unique_sample_id_header = "unique_sample_id"


# query braod list
def map_clinical_phase_to_number(phase):
    match (str(phase)):
        case "" | "None" | "NaN" | "nan" | "np.nan" | "0" | "No Development Reported":
            return 0
        case "Preclinical":
            return 0.5
        case "1" | "1.0" | "Phase 1":
            return 1
        case "2" | "2.0" | "Phase 2":
            return 2
        case "3" | "3.0" | "Phase 3":
            return 3
        case "4" | "4.0" | "Phase 4" | "Launched":
            return 4
        case _:
            return phase


def get_clinical_phase_description(number):
    match (str(number)):
        case "" | "None" | "NaN" | "nan" | "np.nan" | "0" | "0.0":
            return ""
        case "0.5":
            return "Preclinical"
        case "1.0" | "1":
            return "Phase 1"
        case "2.0" | "2":
            return "Phase 2"
        case "3.0" | "3":
            return "Phase 3"
        case "4.0" | "4":
            return "Launched"
        case _:
            return number


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

    ensure_synonyms_column(df)

    # Query pubchem by name and CAS
    if query_pubchem_by_name:
        df = pubchem_search_structure_by_name(df)

    # get mol from smiles or inchi
    # calculate all identifiers from mol - exact_mass, ...
    if calc_identifiers:
        clean_structure_add_mol_id_columns(df)

    # drop duplicates because PubChem name search might generate new rows for conflicting smiles structures
    df = drop_duplicates_by_structure_rowid(df)

    # add new columns for cross references to other databases
    if query_unichem:
        unichem_df = unichem_client.search_all_xrefs(df)
        unichem_client.save_unichem_df(metadata_file, unichem_df)
        df = unichem_client.extract_ids_to_columns(unichem_df, df)

    if query_pubchem_by_structure:
        df = pubchem_search_by_structure(df)

    # extract ids like the UNII, ...
    df = ensure_synonyms_column(df)
    df = extract_synonym_ids(df)

    if query_chembl:
        df = chembl_search_id_and_inchikey(df)

    if query_broad_list:
        df = broad_list_search(df)

    if query_drugbank_list:
        df = drugbank_search_add_columns(df)

    if query_drugcentral:
        df = drugcentral_search(df)

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

    # GNPS cached version
    if query_npclassifier:
        df = apply_np_classifier(df)

    # GNPS cached version
    if query_classyfire:
        df = apply_classyfire(df)

    # export metadata file
    logging.info("Exporting to file %s", out_file)
    if metadata_file.endswith(".tsv"):
        df.to_csv(out_file, sep="\t", index=False)
    else:
        df.to_csv(out_file, sep=",", index=False)


def drop_duplicates_by_structure_rowid(df):
    try:
        id_columns = ["inchikey", "row_id"]
        df = df.drop_duplicates(id_columns, keep="first").sort_index()
    except:
        pass
    return df


if __name__ == "__main__":
    cleanup_file(r"data/library/test_metadata.tsv", lib_id="pluskal_nih",
                 query_pubchem_by_name=True,
                 # need local files
                 query_broad_list=True, query_drugbank_list=True, query_drugcentral=True)
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
