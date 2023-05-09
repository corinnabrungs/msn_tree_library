import pandas as pd
from pandas import DataFrame
from datetime import date
from tqdm import tqdm
import logging
import os

tqdm.pandas()
logging.basicConfig(format="%(asctime)s - %(message)s", level=logging.DEBUG)

base_filename_header = "base_filename"


def main():
    metadata_file = r"data/nih/nih_library_new_headers.tsv"
    data_filepath = r"C:\Xcalibur\data\Corinna_Brungs"

    # define all variables
    inject_volume_mul = 3
    lib_id = "pluskal_nih"
    method_suffix = "TESTMETHOD"  # is added to the end of the path and file names
    instrument_method_positive = r"C:\Xcalibur\methods\Corinna_Brungs\IT_acquisition\IT_100AG_MS5_POS_mz115-2000"
    instrument_method_negative = r"C:\Xcalibur\methods\Corinna_Brungs\IT_acquisition\IT_100AG_MS5_NEG_mz115-2000"

    # plates are inserted into x compartment
    # if plate is not named in the metadata table by plate_id_header column - leave the plate_id empty
    plates_in_autosampler_location = [
        # plate_id, location
        ("15", "B"),
    ]

    create_orbitrap_sequence(metadata_file, data_filepath, instrument_method_positive, instrument_method_negative,
                             lib_id, method_suffix, plates_in_autosampler_location, inject_volume_mul=inject_volume_mul)


def create_orbitrap_sequence(metadata_file, data_filepath: str, instrument_method_positive: str | None,
                             instrument_method_negative: str | None, lib_id: str, method_suffix: str,
                             plates_in_autosampler_location: list, unique_id_header="unique_sample_id",
                             plate_id_header="plate_id", well_header="well_location", inject_volume_mul=3):
    """
    Creates sequences for orbitrap instruments
    :param metadata_file: the metadata file that contains the well_location, plate_id, and unique_sample_id columns
    :param data_filepath: path to store acquired data to
    :param instrument_method_positive: positive mode method - or None to skip
    :param instrument_method_negative: negative mode method - or None to skip
    :param lib_id: defines the compound library
    :param method_suffix: defines the method, e.g., MSn, IT, HCD, ...
    :param plates_in_autosampler_location: list of tuples plate_id, location as tuples ("15", "B"),
    :param unique_id_header: defines a unique sample id. Must not end or start with a number so that contains matches are unique even for A1 and A10
    :param plate_id_header: plate id column. Plate ids can be any string or number
    :param well_header: defines the column with well locations, e.g., A1
    :param inject_volume_mul: micro liter injection volume
    :return:
    """
    current_date = date.today().strftime("%Y%m%d")
    # NO NEED TO CHANGE ANYTHING BELOW
    # final values
    data_filepath = os.path.join(data_filepath, lib_id, f"{current_date}_{method_suffix}")
    dataframes = []
    metadata_df = load_metadata_df(metadata_file, well_header, unique_id_header, method_suffix)
    for plate_id, plate_location in plates_in_autosampler_location:
        plate_df = filter_metadata_by_plate_id(metadata_df, plate_id, plate_id_header)

        sequence_file = f"data/Sequence/{current_date}_seq_rack_{plate_location}_{lib_id}_{plate_id}_{method_suffix}"

        df = _create_orbitrap_sequence(plate_df, sequence_file, data_filepath, well_header, plate_location,
                                       instrument_method_positive, instrument_method_negative, inject_volume_mul)
        dataframes.append(df)
    concat = pd.concat(dataframes)
    plates_str = "_".join(["{}in{}".format(plate_id, loc) for plate_id, loc in plates_in_autosampler_location])
    sequence_file = f"data/Sequence/{current_date}_{plates_str}_seq_combined.csv"
    write_thermo_sequence(sequence_file, concat)


def _create_orbitrap_sequence(metadata_df: DataFrame, sequence_file, data_filepath, well_header, plate_location,
                              instrument_method_positive=None, instrument_method_negative=None,
                              inject_volume_mul=3) -> DataFrame:
    """
    Creates Orbitrap sequence for positive and negative mode

    :param data_filepath: data acquisition file path. polarity will be added
    :param sequence_file: the base sequence file to export. polarity and file type csv will be added automatically
    :param metadata_df: a dataframe that is already filtered to only contain samples from a single plate
    :param well_header: getting the well number of the final plate, e.g., A1
    :param instrument_method_positive: instrument method in positive mode
    :param instrument_method_negative: instrument method in negative mode
    :param plate_location: position in the autosampler
    Sequence for Plate 2 and Position Green
    :param inject_volume_mul: injection volume in micro liter
    :return:
    """
    if not instrument_method_positive and not instrument_method_negative:
        raise ValueError("Provide at least one method file for positive or negative")

    instrument_methods = [instrument_method_positive, instrument_method_negative]
    polarities = ["positive", "negative"]
    dataframes = []
    for polarity, instrument_method in zip(polarities, instrument_methods):
        if not instrument_method:
            continue

        seq_df = DataFrame()
        seq_df["File Name"] = [f"{base_filename}_{polarity}" for base_filename in metadata_df[base_filename_header]]
        #
        seq_df["Path"] = f"{data_filepath}_{polarity}"
        seq_df["Instrument Method"] = instrument_method
        seq_df["Position"] = ["{}:{}".format(plate_location, well) for well in metadata_df[well_header]]

        seq_df["Inj Vol"] = inject_volume_mul
        seq_df["Dil Factor"] = 1
        seq_df = seq_df.drop_duplicates()

        csv_file = f"{sequence_file}_{polarity}.csv"
        write_thermo_sequence(csv_file, seq_df)
        dataframes.append(seq_df)
    return pd.concat(dataframes)


def load_metadata_df(metadata_file, well_header, unique_id_header, method_suffix) -> DataFrame:
    logging.info("Will run on %s", metadata_file)
    # import df
    if metadata_file.endswith(".tsv"):
        df = pd.read_csv(metadata_file, sep="\t")
    else:
        df = pd.read_csv(metadata_file, sep=",")

    if well_header not in df.columns:
        raise ValueError(f"No column named {well_header} with the well number, e.g., A1")
    if unique_id_header not in df.columns:
        raise ValueError(
            f"No column named {unique_id_header} with unique sample ids. Run metadata clean up that generates a unique id, e.g., lib_plate1_A1_id (note the _id at the end and the prefix that make sure that wells like A1 do not match to A10)")

    current_date = date.today().strftime("%Y%m%d")
    df[base_filename_header] = ["{}_{}_{}".format(current_date, unique_id, method_suffix) for unique_id in
                                df[unique_id_header]]
    return df


def filter_metadata_by_plate_id(metadata_df: DataFrame, plate_id: str, plate_id_header: str) -> DataFrame:
    """
    Filter metadata_df to only contain entries for the current plate_id
    DONT apply filter if plate_id is empty - instead use full metadata_df
    :param metadata_df:
    :param plate_id:
    :param plate_id_header:
    :return:
    """
    if plate_id:
        if plate_id_header not in metadata_df.columns:
            raise ValueError(
                f"Plate id filter was {plate_id} but there was no column in the metadata file for the plates. Name: {plate_id_header}")
        else:
            plate_df = metadata_df[metadata_df[plate_id_header].astype(str) == plate_id]
    else:
        plate_df = metadata_df
    return plate_df


def write_thermo_sequence(csv_file, df: DataFrame):
    df.to_csv(csv_file, index=False)
    # df.to_csv("data/nih/{}_uniqueID.csv".format(lib_id), index=False)
    # Adding the first line as needed by Xcalibur sequence
    with open(csv_file, 'r') as original:
        data = original.read()
    with open(csv_file, 'w') as modified:
        modified.write("Bracket Type=4,\n" + data)

    logging.info(f"Saved new sequence to: {csv_file}")


if __name__ == "__main__":
    try:
        main()
    except:
        logging.exception("Could not create sequence")
        exit(1)
    exit(0)
