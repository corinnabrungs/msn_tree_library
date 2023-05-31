import logging
from pathlib import Path

import pandas as pd
import numpy as np
from pandas._typing import IndexLabel


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
    elif file.endswith(('.parquet', '.parquet.gz', '.parquet.gzip')):
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
    elif out_file.endswith('.csv'):
        df.to_csv(out_file, sep=",", index=False)
    elif out_file.endswith('.parquet'):
        df.to_parquet(out_file)
    elif out_file.endswith('.parquet.gz') or out_file.endswith('.parquet.gzip'):
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


def add_column_prefix(df: pd.DataFrame, prefix: str, columns_to_rename=None, columns_to_keep=None) -> pd.DataFrame:
    """
    Add prefix to all columns in "to rename" but never to those in "to keep"
    :param df: dataframe
    :param prefix: prefix+old column name
    :return: the input df
    """
    if columns_to_rename is None:
        columns_to_rename = df.columns
    if isinstance(columns_to_rename, str):
        columns_to_rename = [columns_to_rename]

    if columns_to_keep is None:
        columns_to_keep = []
    if isinstance(columns_to_keep, str):
        columns_to_keep = [columns_to_keep]

    df.columns = [col if col in columns_to_keep else prefix + col for col in columns_to_rename]
    return df


def remove_empty_strings(df: pd.DataFrame, columns) -> pd.DataFrame:
    if isinstance(columns, str):
        columns = [columns]

    for col in columns:
        if col in df.columns:
            df[col] = [None if isinstance(v, str) and len(v) == 0 else v for v in df[col]]

    return df


def remove_empty_lists_values(df: pd.DataFrame, columns) -> pd.DataFrame:
    if isinstance(columns, str):
        columns = [columns]

    for col in columns:
        if col in df.columns:
            df[col] = [None if isinstance(v, list) and len(v) == 0 else v for v in df[col]]

    return df


def join_unique(values, separator: str = "; "):
    try:
        if isinstance(values, str):
            return values
        return separator.join(get_unique_list(values, False))
    except:
        return None


def groupby_join_unique_values(df: pd.DataFrame, columns, as_lists: bool = False,
                               str_separator: str = "; ") -> pd.DataFrame:
    if isinstance(columns, str):
        columns = [columns]

    if as_lists:
        return df.groupby(columns).agg(lambda col: list(set(col))).reset_index()
    else:
        return df.groupby(columns).agg(lambda col: join_unique(col, str_separator)).reset_index()


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


def left_merge_retain_index(main_index_df: pd.DataFrame, other_df: pd.DataFrame,
                            on: IndexLabel | None = None) -> pd.DataFrame:
    """
    Merge other_df into main_df with left join and keep main_df.index
    :param main_index_df: provides index
    :param other_df: is merged into main
    :param on: columns, or if none use index
    :return: merged frame
    """
    return main_index_df.merge(other_df, how="left", sort=False, on=on).set_index(main_index_df.index)


def get_first_value_or_else(df: pd.DataFrame, column: str, default=None):
    return next((v for v in df[column]), default)


def get_or_else(row, key, default=None):
    return row[key] if key in row and notnull(row[key]) else default


def get_unique_list(input_list, ignore_case: bool = True):
    seen = set()
    if ignore_case:
        return [x for x in input_list if notnull(x) and x.lower() not in seen and not seen.add(x.lower())]
    else:
        return [x for x in input_list if notnull(x) and x not in seen and not seen.add(x)]


def get_unique_dict(df, column):
    """
    Dict with unique keys and no values
    :param df: input data frame
    :return: dict(key, None)
    """
    return dict.fromkeys(df[column])


def create_missing_columns(df: pd.DataFrame, required_cols):
    if isinstance(required_cols, str):
        required_cols = [required_cols]
    for col in required_cols:
        if col not in df.columns:
            df[col] = None
    return df


def make_str_floor_to_int_number(df: pd.DataFrame, columns) -> pd.DataFrame:
    """
    Works in place. changes all input columns to str, splits at . and takes the first part to ensure integers.
    "1.0" will be "1"
    :param df: data frame input
    :param columns: str or list of columns
    :return: the input data frame
    """
    if isinstance(columns, str):
        columns = [columns]
    for col in columns:
        df[col] = df[col].astype("str").str.split('.').str[0]
    return df


def divide_chunks(items, chunk_size):
    """

    :return: list of chunks of specified size
    """
    return [items[i:i + chunk_size] for i in range(0, len(items), chunk_size)]
