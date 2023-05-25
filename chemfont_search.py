import metadata_cleanup as cleanup
import logging
import pandas as pd
from tqdm import tqdm

import chemfont_postgresql_query as chemfont_query
import pandas_utils

tqdm.pandas()
logging.basicConfig(format="%(asctime)s - %(message)s", level=logging.DEBUG)


def chembl_availability(df) -> pd.DataFrame:
    compounds = [cleanup.client.get_chembl_mol(chembl_id, inchi_key) for chembl_id, inchi_key in
                 tqdm(zip(df["chembl_id"], df["inchi_key"]))]
    df["availability"] = [compound["availability_type"] if notnull(compound) else np.NAN for compound in compounds]
    return df


def chemfont_search(df):
    prefix = "chemfont_"
    try:
        chemfont_query.connect()
        logging.info("Searching in DrugCentral")
        results = df.progress_apply(lambda row: chemfont_query.chemfont_for_row(row), axis=1)
        # results = [drugcentral_query.drugcentral_for_row(row) for _,row in df.iterrows()]
        # results = [drugcentral_query.drugcentral_postgresql(inchikey, split_inchikey) for inchikey, split_inchikey in
        #          tqdm(zip(df["inchi_key"], df["split_inchi_key"]))]
        logging.info("DrugCentral search done")

        # row[1] is the data row[0] is the columns
        first_entry_columns = next((row[0] for row in results if notnull(row[0])), [])
        columns = [col.name for col in first_entry_columns]
        elements = len(columns)
        data = [row[1] if notnull(row[1]) else (None,) * elements for row in results]
        chemfont_df = pd.DataFrame(data=data, columns=columns, index=df.index)
        chemfont_df = chemfont_df.add_prefix(prefix)
        return pd.concat([df, chemfont_df], axis=1)

    finally:
        chemfont_query.deconnect()


def chemfont_file(metadata_file):
    logging.info("Will run on %s", metadata_file)
    out_file = pandas_utils.add_filename_suffix(metadata_file, "chemfont")

    if metadata_file.endswith(".tsv"):
        df = pd.read_csv(metadata_file, sep="\t")
    else:
        df = pd.read_csv(metadata_file, sep=",")

    logging.info("Search ChEMBL by chemblid or inchikey")
    df = chembl_availability(df)

    logging.info("Search chemfont by external identifier or inchikey")
    df = chemfont_search(df)

    logging.info("Exporting to file %s", out_file)
    if metadata_file.endswith(".tsv"):
        df.to_csv(out_file, sep="\t", index=False)
    else:
        df.to_csv(out_file, sep=",", index=False)


if __name__ == "__main__":
    chemfont_file(r"data\test_metadata.tsv")
