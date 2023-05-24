from pathlib import Path

import pandas as pd


def isnull(o):
    return str(o) == '<NA>' or o == "N/A" or o == "NA" or o == None or pd.isnull(o)


def notnull(o):
    return not isnull(o)


def read_dataframe(file):
    if file.endswith(".tsv"):
        df = pd.read_csv(file, sep="\t")
    elif file.endswith('.csv'):
        df = pd.read_csv(file)
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


def add_filename_suffix(filename, suffix, format_override=None):
    """

    :param filename: original file
    :param suffix: is added to the filename
    :param format_override: None will use the original file suffix and otherwise specify changed suffix
    :return: filename_suffix.format
    """
    p = Path(filename)
    file_format = p.suffix if format_override is None else format_override
    return "{0}_{1}{2}".format(Path.joinpath(p.parent, p.stem), suffix, file_format)


def combine_dfs_fill_missing_values(target: pd.DataFrame, source: pd.DataFrame) -> pd.DataFrame:
    """
    Only missing values are replaced. Column names need to match between the target and source
    :param target: has priority. only missing values are filled
    :param source: source to fill values
    :return: the filled dataframe
    """
    return target.combine_first(source)  # alternative df.combine_first


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
