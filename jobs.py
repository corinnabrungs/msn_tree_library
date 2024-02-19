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
    (
        r"C:\git\msn_library\data\library\mce_new_run_parallel\mce_library_all.tsv",
        "mce",
    ),
    # (r"C:\git\msn_library\data\nina\one_hundred_drug_test.csv", "nina_drug"),
    # (r"C:\git\msn_library\data\nina\one_hundred_drug_test.csv", "nina_drug"),
    ###
    # (r"C:\git\msn_library\data\nina\reframe\nina_reframe_chunk0.csv", "nina_reframe"),
    # (r"C:\git\msn_library\data\nina\reframe\nina_reframe_chunk1.csv", "nina_reframe"),
    # (r"C:\git\msn_library\data\nina\reframe\nina_reframe_chunk2.csv", "nina_reframe"),
    # (r"C:\git\msn_library\data\nina\reframe\nina_reframe_chunk3.csv", "nina_reframe"),
    # (r"C:\git\msn_library\data\nina\reframe\nina_reframe_chunk4.csv", "nina_reframe"),
    # (r"C:\git\msn_library\data\nina\reframe\nina_reframe_chunk5.csv", "nina_reframe"),
    # (r"C:\git\msn_library\data\nina\reframe\nina_reframe_chunk6.csv", "nina_reframe"),
    # (r"C:\git\msn_library\data\nina\reframe\nina_reframe_chunk7.csv", "nina_reframe"),
    # (r"C:\git\msn_library\data\nina\reframe\nina_reframe_chunk8.csv", "nina_reframe"),
    # (r"C:\git\msn_library\data\nina\reframe\nina_reframe_chunk9.csv", "nina_reframe"),
]

if __name__ == "__main__":
    # make sure prefect server runs and its locally deployed
    # keep metadata cleanup prefect running

    cfg = MetadataCleanupConfig(
        query_npatlas=False,
        query_broad_list=True,
        query_drugbank_list=True,
        query_drugcentral=True,
        query_lotus=True,
    )

    for metadata_file, lib_id in files_and_lib_ids:
        try:
            run_async(
                cfg=cfg,
                metadata_file=metadata_file,
                lib_id=lib_id,
                flow_method=cleanup_file_chunked,
                use_cached_parquet_file=True,
                full_iterations=2,
            )
        except:
            logging.exception("Exception in flow")

    exit(0)
