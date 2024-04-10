import logging
import os
import sys

import prefect
from prefect.deployments import run_deployment

import metadata_cleanup_prefect

from concurrent.futures import ThreadPoolExecutor
from concurrent.futures import wait

from metadata_cleanup_prefect import (
    cleanup_file,
    MetadataCleanupConfig,
    cleanup_file_chunked,
    run_async,
)

logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(__name__)

files_and_lib_ids = [
    (r"examples\test_metadata.tsv", "test"),
    # (r"examples\test_metadata_medium.tsv", "test_medium"),
    # (r"examples\test_metadata_small.tsv", "test_small"),
]

if __name__ == "__main__":
    # make sure prefect server runs and its locally deployed
    # keep metadata cleanup prefect running

    cfg = MetadataCleanupConfig(
        query_npatlas=True,
        query_broad_list=False,
        query_drugbank_list=False,
        query_drugcentral=False,
        query_lotus=False,
        query_classyfire=True,
        query_dictionary_np=False,
    )

    for metadata_file, lib_id in files_and_lib_ids:
        try:
            # run_async(
            #     cfg=cfg,
            #     metadata_file=metadata_file,
            #     lib_id=lib_id,
            #     flow_method=cleanup_file_chunked,
            #     use_cached_parquet_file=True,
            #     full_iterations=1,
            # )
            run_async(
                cfg=cfg,
                metadata_file=metadata_file,
                lib_id=lib_id,
                flow_method=cleanup_file,
                use_cached_parquet_file=True,
                full_iterations=1,
            )
        except:
            logging.exception("Exception in flow")

    exit(0)
