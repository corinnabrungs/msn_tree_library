import dataclasses
import logging
from dataclasses import dataclass
from pathlib import Path

import prefect
from prefect.deployments import run_deployment
from tqdm import tqdm
import lotus_client
from dictionary_of_np_client import dictionary_of_np_search
import pandas_utils as pu
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
from pandas_utils import (
    update_dataframes,
    create_missing_columns,
    read_dataframe,
    check_if_chunks_available,
    add_filename_suffix,
)
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


@dataclass
class MetadataCleanupConfig:
    plate_id_header: str = "plate_id"
    well_header: str = "well_location"
    query_pubchem_by_cid: bool = True
    query_pubchem_by_name: bool = True
    calc_identifiers: bool = True
    query_unichem: bool = True
    query_pubchem_by_structure: bool = True
    query_chembl: bool = True
    query_npclassifier: bool = True
    query_classyfire: bool = True
    query_npatlas: bool = False
    query_broad_list: bool = False
    query_drugbank_list: bool = False
    query_drugcentral: bool = False
    query_lotus: bool = False
    query_dictionary_np: bool = False


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


@task(name="dictionary_np_search")
def dictionary_of_np_search_prefect(df):
    return dictionary_of_np_search(df)


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
    name="Metadata cleanup-chunked",
    version="0.2.0",
    flow_run_name="{lib_id}:{metadata_file}",
)
def cleanup_file_chunked(
    cfg: MetadataCleanupConfig,
    metadata_file: str,
    lib_id: str,
    use_cached_parquet_file: bool = True,
    n_thread=1,
    min_chunk_size=250,
    max_chunk_size=1000000,
    full_iterations=1,
):
    if not use_cached_parquet_file:
        pu.delete_chunks(metadata_file)
    # check if already chunks
    has_chunks = check_if_chunks_available(metadata_file)
    if not has_chunks:
        # read df
        df = read_dataframe(metadata_file)
        # >1000 split into chunks of 1000 rows in .parquet
        if len(df) > max_chunk_size or n_thread > 1:
            chunks = pu.save_chunks(
                df, metadata_file, n_thread, min_chunk_size, max_chunk_size
            )
            has_chunks = len(chunks) > 1  # if 1 then original df

    if not has_chunks:
        # run one job
        run_async(
            cfg=cfg,
            metadata_file=metadata_file,
            lib_id=lib_id,
            flow_method=cleanup_file,
            use_cached_parquet_file=use_cached_parquet_file,
            full_iterations=full_iterations,
        )
    else:
        counter = 0
        while True:  # loop chunks
            try:
                file = pu.add_filename_suffix(
                    metadata_file, f"chunk{counter}", ".parquet"
                )
                if not Path(file).is_file():
                    break
                run_async(
                    cfg=cfg,
                    metadata_file=file,
                    lib_id=lib_id,
                    flow_method=cleanup_file,
                    use_cached_parquet_file=use_cached_parquet_file,
                    full_iterations=full_iterations,
                )
                counter += 1
            except:
                break


def run_async(
    cfg: MetadataCleanupConfig,
    metadata_file: str,
    lib_id: str,
    flow_method,
    deployment_name="local-deploy",
    use_cached_parquet_file: bool = True,
    current_iteration: int | None = None,
    full_iterations: int = 2,
):
    params = {
        "cfg": dataclasses.asdict(cfg),
        "metadata_file": metadata_file,
        "lib_id": lib_id,
        "use_cached_parquet_file": use_cached_parquet_file,
        "full_iterations": full_iterations,
    }
    if current_iteration is not None:
        params["current_iteration"] = current_iteration

    run_deployment(
        f"{flow_method.name}/{deployment_name}",
        parameters=params,
        timeout=0,  # non blocking, to not wait for flow to finish
    )


@flow(
    name="Metadata cleanup", version="0.2.0", flow_run_name="{lib_id}:{metadata_file}"
)
def cleanup_file(
    cfg: MetadataCleanupConfig,
    metadata_file: str,
    lib_id: str,
    use_cached_parquet_file: bool = True,
    current_iteration: int = 1,
    full_iterations: int = 1,
):
    if current_iteration > 1:
        use_cached_parquet_file = True

    logging.info("Will run on %s", metadata_file)
    df = extract_prepare_input_data_prefect(
        metadata_file,
        lib_id,
        cfg.plate_id_header,
        cfg.well_header,
        use_cached_parquet_file,
    )

    if cfg.query_pubchem_by_cid:
        df = pubchem_search_parent_by_cid_prefect(df, apply_structures=True)

    # Query pubchem by name and CAS
    if cfg.query_pubchem_by_name:
        df = pubchem_search_structure_by_name_prefect(df)

    save_intermediate_parquet_prefect(df, metadata_file)

    # get mol from smiles or inchi
    # calculate all identifiers from mol - monoisotopic_mass, ...
    if cfg.calc_identifiers:
        df = clean_structure_add_mol_id_columns_prefect(df)

    if cfg.query_pubchem_by_structure:
        df = pubchem_search_by_structure_prefect(df)
        save_intermediate_parquet_prefect(df, metadata_file)

    if cfg.query_pubchem_by_cid:
        df = pubchem_search_parent_by_cid_prefect(df, apply_structures=True)
        if cfg.calc_identifiers:
            df = clean_structure_add_mol_id_columns_prefect(df)
        save_intermediate_parquet_prefect(df, metadata_file)

    # drop duplicates because PubChem name search might generate new rows for conflicting smiles structures
    df = drop_duplicates_by_structure_rowid_reset_index(df)
    save_intermediate_parquet_prefect(df, metadata_file)

    # structures are now fetched. Run things in parallel
    # run in parallel
    tasks = []

    # GNPS cached version
    if cfg.query_npclassifier:
        tasks.append(apply_np_classifier_prefect.submit(df))

    # GNPS cached version
    if cfg.query_classyfire:
        tasks.append(apply_classyfire_prefect.submit(df))

    if cfg.query_npatlas:
        tasks.append(search_np_atlas_prefect.submit(df))

    if cfg.query_lotus:
        tasks.append(search_lotus_prefect.submit(df))

    if cfg.query_dictionary_np:
        df = dictionary_of_np_search_prefect(df)

    # add new columns for cross references to other databases
    if cfg.query_unichem:
        # xrefs are needed for other steps so run sequential here
        df = search_all_unichem_xrefs_prefect(df, metadata_file)
        save_intermediate_parquet_prefect(df, metadata_file)

    # extract ids like the UNII, ...
    df = ensure_synonyms_column(df)
    df = extract_synonym_ids(df)

    if cfg.query_chembl:
        df = chembl_search_id_and_inchikey_prefect(df)
        save_intermediate_parquet_prefect(df, metadata_file)

    if cfg.query_broad_list:
        df = broad_list_search_prefect(df)

    if cfg.query_drugbank_list:
        df = drugbank_search_add_columns_prefect(df)

    if cfg.query_drugcentral:
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

    if current_iteration < full_iterations:
        # run task again to fill all values and retry failed services
        # run asynch so that this flow can finish and give way to the new call
        run_async(
            cfg=cfg,
            metadata_file=metadata_file,
            lib_id=lib_id,
            flow_method=cleanup_file,
            use_cached_parquet_file=True,
            current_iteration=current_iteration + 1,
            full_iterations=full_iterations,
        )


def full_cleanup_file_chunked(
    metadata_file, lib_id, use_cached_parquet_file: bool = True
):
    try:
        cfg = MetadataCleanupConfig(
            query_npatlas=False,
            query_broad_list=True,
            query_drugbank_list=True,
            query_drugcentral=True,
            query_lotus=True,
        )
        run_async(
            cfg=cfg,
            metadata_file=metadata_file,
            lib_id=lib_id,
            flow_method=cleanup_file_chunked,
            use_cached_parquet_file=use_cached_parquet_file,
            full_iterations=2,
        )
    except:
        logging.exception("Exception in flow")
        # exit(1)


if __name__ == "__main__":
    prefect.serve(
        cleanup_file_chunked.to_deployment(
            "local-deploy",
            work_pool_name="local-work",
        ),
        cleanup_file.to_deployment(
            "local-deploy",
            work_pool_name="local-work",
        ),
    )
