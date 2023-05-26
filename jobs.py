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
    # (r"data\nih\nih_library_test.csv", "pluskal_nih"),
    (r"data\nih\nih_library_new_headers.tsv", "pluskal_nih"),
    # (r"examples\test_metadata.tsv", "test"),
    # (r"examples\test_metadata_small.tsv", "test"),
    # (r"data\library\mce_library.tsv", "pluskal_mce"),
    # (r"data\library\mce_library_add_compounds.tsv", "pluskal_mce"),
]

if __name__ == "__main__":
    import subprocess

    for input_file, lib_id in files_and_lib_ids:
        metadata_cleanup_prefect.full_cleanup_file(input_file, lib_id, flow_run_name="pluskal_nih")

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
