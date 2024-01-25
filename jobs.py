import logging
import os
import sys

import prefect

import metadata_cleanup_prefect

from concurrent.futures import ThreadPoolExecutor
from concurrent.futures import wait

logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(__name__)

files_and_lib_ids = [
    # (r"examples\test_metadata.tsv", "test"),
    (r"examples\test_metadata_small.tsv", "test"),
]

if __name__ == "__main__":
    import subprocess

    for input_file, lib_id in files_and_lib_ids:
        # use_cached_parquet_file = False cleans up the parquet files - use True instead
        try:
            metadata_cleanup_prefect.cleanup_file(
                input_file,
                lib_id,
                use_cached_parquet_file=True,
                query_pubchem_by_name=True,
                query_pubchem_by_cid=True,
                query_unichem=True,
                query_chembl=True,
                query_npatlas=True,
                query_classyfire=True,
                # need local files
                query_broad_list=False,
                query_drugbank_list=False,
                query_drugcentral=False,
                query_lotus=False,
            )
        except:
            logging.exception("Exception in flow")
            # exit(1)

        # metadata_cleanup_prefect.full_cleanup_file(input_file, lib_id, use_cached_parquet_file=True)
        # metadata_cleanup_prefect.add_lotus_flow(input_file, lib_id)

        # cmd = f"python metadata_cleanup_prefect.py {input_file} --lib_id {lib_id}"
        # logging.info("Starting file "+input_file)
        # os.system(cmd)
        # subprocess.call(cmd, shell=True)

        # p = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
        # out, err = p.communicate()
        # result = out.split('\n')
        # for lin in result:
        #     if not lin.startswith('#'):
        #         print(lin)

    # parallel_queries = 8
    # with ThreadPoolExecutor(parallel_queries) as executor:
    #     futures = [
    #         executor.submit(
    #             metadata_cleanup_prefect.full_cleanup_file,
    #             input_file,
    #             lib_id
    #         )
    #         for input_file, lib_id in files_and_lib_ids
    #     ]
    #
    #     wait(futures)
    # exit with OK
    sys.exit(0)
