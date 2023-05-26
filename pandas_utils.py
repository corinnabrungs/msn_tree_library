import logging
from pathlib import Path

import pandas as pd
import numpy as np


def isnull(o):
    try:
        if o is None:
            return True

        so = str(o).lower()
        if so == '<na>' or so == "n/a" or so == "na" or so == "nan" or so == "nat":
            return True

        result = pd.isnull(o)
        if isinstance(result, bool):
            return result

        # otherwise we are looking at a series or list which is not None
        return False
    except:
        return False


def notnull(o):
    return not isnull(o)


def read_dataframe(file):
    if file.endswith(".tsv"):
        df = pd.read_csv(file, sep="\t")
    elif file.endswith('.csv'):
        df = pd.read_csv(file)
    elif file.endswith('.parquet', '.parquet.gz', '.parquet.gzip'):
        df = pd.read_parquet(file)
    elif file.endswith('.feather'):
        df = pd.read_feather(file)
    elif file.endswith('.xls', 'xlsx'):
        df = pd.read_excel(file)
    elif file.endswith('.xml'):
        df = pd.read_xml(file)
    elif file.endswith('.json'):
        df = pd.read_json(file)
    elif file.endswith('.sql'):
        df = pd.read_sql(file)
    elif file.endswith('.hdf'):
        df = pd.read_hdf(file)
    else:
        raise ValueError(f'Unsupported filetype: {file}')
    return df


def save_dataframe(df, out_file):
    logging.info("Exporting to file %s", out_file)
    if out_file.endswith(".tsv"):
        df.to_csv(out_file, sep="\t", index=False)
    elif out_file.endswith('.parquet'):
        df.to_parquet(out_file)
    elif out_file.endswith('.parquet.gz', '.parquet.gzip'):
        df.to_parquet(out_file, compression='gzip')
    elif out_file.endswith('.feather'):
        df.to_feather(out_file)
    else:
        df.to_csv(out_file, sep=",", index=False)


def get_parquet_file(metadata_file, gzip=False):
    if gzip:
        return replace_format(metadata_file, ".parquet.gzip")
    return replace_format(metadata_file, ".parquet")


def replace_format(filename, format_override):
    """

    :param filename: original file
    :param format_override: change file format
    :return: filename.format
    """
    format_override = format_override if format_override.startswith(".") else "." + format_override
    p = Path(filename)
    return "{0}{1}".format(Path.joinpath(p.parent, p.stem), format_override)


def add_filename_suffix(filename, suffix, format_override=None):
    """

    :param filename: original file
    :param suffix: is added to the filename
    :param format_override: None will use the original file suffix and otherwise specify changed suffix
    :return: filename_suffix.format
    """
    p = Path(filename)
    file_format = p.suffix if format_override is None else format_override
    file_format = file_format if file_format.startswith(".") else "." + file_format
    return "{0}_{1}{2}".format(Path.joinpath(p.parent, p.stem), suffix, file_format)


def remove_empy_strings(df: pd.DataFrame, columns) -> pd.DataFrame:
    if isinstance(columns, str):
        columns = [columns]

    for col in columns:
        if col in df.columns:
            df[col] = [v if isinstance(v, str) and len(v) > 0 else None for v in df[col]]

    return df


def update_dataframes(newdf: pd.DataFrame, olddf: pd.DataFrame) -> pd.DataFrame:
    """
    Merges old df into newdf so that new columns are added and filled. NA values in newdf are filled with olddf values.
    :param newdf: new dataframe, can be a filtered slice of the old df
    :param olddf: the old original dataframe
    :return: new dataframe
    """
    return combine_dfs_fill_missing_values(newdf, olddf)


def combine_dfs_fill_missing_values(target: pd.DataFrame, source: pd.DataFrame) -> pd.DataFrame:
    """
    Only missing values are replaced. Column names need to match between the target and source
    :param target: has priority. only missing values are filled
    :param source: source to fill values
    :return: the filled dataframe
    """
    return target.combine_first(source)  # alternative df.combine_first


def get_first_value_or_else(df: pd.DataFrame, column: str, default=None):
    return next((v for v in df[column]), default)


def get_or_else(row, key, default=None):
    return row[key] if key in row and notnull(row[key]) else default


def get_unique_list(input_list):
    seen = set()
    unique = [x for x in input_list if x.lower() not in seen and not seen.add(x.lower())]
    return unique


def get_unique_dict(df, column):
    """
    Dict with unique keys and no values
    :param df: input data frame
    :return: dict(key, None)
    """
    return dict.fromkeys(df[column])
