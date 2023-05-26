import datetime as dt

import pandas as pd

from pandas_utils import isnull


def parse_str(iso_datetime):
    return dt.datetime.fromisoformat(iso_datetime)


def iso_datetime_now() -> str:
    """
    Example 2023-05-26T11:45:54+02:00
    :return: day time timezone as yyyy-mm-ddThh:mm:ss+hh:mm
    """
    return datetime_now().isoformat()


def datetime_now() -> dt.datetime:
    """
    :return: time-zoned datetime with second precision
    """
    return dt.datetime.now().replace(microsecond=0).astimezone()


def check_within_timedelta_from_now(old_datetime: dt.datetime | str | None, deltatime: dt.timedelta) -> bool:
    """
    Checks if old datetime is within now - deltatime
    :param old_datetime: the checked old time or None
    :param deltatime: a timedelta
    :return: True if (now - old_datetime) <= deltatime; False if old_datetime is None
    """
    if isnull(old_datetime):
        return False
    if isinstance(old_datetime, str):
        old_datetime = parse_str(old_datetime)
    return (datetime_now() - old_datetime) <= deltatime


def create_expired_entries_dataframe(df: pd.DataFrame, date_column, refresh_expired_entries_after):
    """
    Check last update datetime and create filtered dataframe of all expired entries.
    Later use pandas_utils.combine_dfs_fill_missing_values
    :param df:
    :param date_column:
    :param refresh_expired_entries_after:
    :return:
    """
    if date_column not in df.columns:
        df[date_column] = None
    dcol = pd.to_datetime(df[date_column])
    expired_entries = [not check_within_timedelta_from_now(date, refresh_expired_entries_after) for date in dcol]
    return df[expired_entries].copy()


def set_date_if_notnull(df: pd.DataFrame, result_column, date_column):
    """
    Set the date_column to the current iso datetime if result_column is notnull
    Changes df in place
    """
    df.loc[df[result_column].notnull(), date_column] = iso_datetime_now()
