import pandas as pd


def combine_polarity(old, new):
    if old == "both":
        return "both"
    match new:
        case "both":
            return new
        case "positive":
            return "both" if old == "negative" else "positive"
        case "negative":
            return "both" if old == "positive" else "negative"
        case _:
            return old


def drop_unique_inchikey_polarity(
    df: pd.DataFrame, results_column: str = "polarity", inchikey_column="inchikey"
) -> pd.DataFrame:
    """
    Merge polarity for all the same inchikeys and then drop duplicates
    """
    unique_dict = {}
    for inchikey, polarity in zip(df[inchikey_column], df["polarity"]):
        oldpolarity = unique_dict.get(inchikey, "missing")
        unique_dict[inchikey] = combine_polarity(oldpolarity, polarity)

    df[results_column] = [unique_dict.get(inchikey) for inchikey in df[inchikey_column]]
    df = (
        df.sort_values(by=["detected"])
        .drop_duplicates([results_column, inchikey_column])
        .sort_index()
    )

    # df[df["inchikey"].duplicated(keep=False)][["inchikey", "polarity", "new_polarity"]]
    return df
