import logging
import os
import sys

import prefect

import metadata_cleanup_prefect

from concurrent.futures import ThreadPoolExecutor
from concurrent.futures import wait
import pandas_utils as pu

logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(__name__)

files_and_num_chunks = [
    (r"C:\git\msn_library\data\nina\reframe\nina_reframe.csv", 10),
]

if __name__ == "__main__":
    for file, nchunks in files_and_num_chunks:
        df = pu.read_dataframe(file)
        pu.save_chunks(df, file, nchunks)
    sys.exit(0)
