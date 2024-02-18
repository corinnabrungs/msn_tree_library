import logging
import os
import sys

import prefect

import metadata_cleanup_prefect

from concurrent.futures import ThreadPoolExecutor
from concurrent.futures import wait

import pandas_utils as pu
from pandas_utils import replace_format, add_filename_suffix

logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(__name__)

files = [
    r"C:\git\msn_library\data\nina\reframe\nina_reframe.csv",
]


if __name__ == "__main__":
    for file in files:
        pu.combine_and_save_chunks(file)
    sys.exit(0)
