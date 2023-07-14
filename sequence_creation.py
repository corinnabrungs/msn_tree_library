import pandas as pd
from pandas import DataFrame
from datetime import date
from tqdm import tqdm
import logging
import os
from dataclasses import dataclass
import pandas_utils as pu

tqdm.pandas()
logging.basicConfig(format="%(asctime)s - %(message)s", level=logging.DEBUG)

base_filename_header = "base_filename"


@dataclass
class InstrumentMethod:
    identifier: str
    path: str


def main():
    metadata_file = r"C:\git\msn_library\data\library\mce_library_all_cleaned.tsv"
    data_filepath = r"D:\Corinna_Brungs"  # storage path for the acquired data

    # define all variables
    inject_volume_mul = 1
    lib_id = "pluskal_mce"
    method_suffix = "rt_ms2"  # is added to the end of the path and file names
    instrument_methods = [
        InstrumentMethod("positive",
                         r"C:\Xcalibur\methods\Corinna_Brungs\RT_prediction\20230615_JLW_method_0_5min_equ_positive_MS2_pos_rcc"),
        # InstrumentMethod("negative", r"Test"),
        # InstrumentMethod("polarity_switching", r"Test"),
    ]

    # plates are inserted into x compartment
    # if plate is not named in the metadata table by plate_id_header column - leave the plate_id empty
    plates_in_autosampler_location = [
        # plate_id, location
        ("1D2", "B"),
        ("1D3", "G"),
    ]

    create_orbitrap_sequence(
        metadata_file, data_filepath, instrument_methods,
        lib_id, method_suffix, plates_in_autosampler_location, inject_volume_mul=inject_volume_mul,
        blank_well="F2", qc_well="F1", blank_qc_autosampler_location="R", blank_every_n_samples=20
    )


def create_orbitrap_sequence(metadata_file, data_filepath: str, instrument_methods: list[InstrumentMethod], lib_id: str,
                             method_suffix: str,
                             plates_in_autosampler_location: list, unique_id_header="unique_sample_id",
                             plate_id_header="plate_id", well_header="well_location", inject_volume_mul=3,
                             blank_well=None, qc_well=None, blank_qc_autosampler_location=None,
                             blank_every_n_samples=20):
    """
    Creates sequences for orbitrap instruments
    :param metadata_file: the metadata file that contains the well_location, plate_id, and unique_sample_id columns
    :param data_filepath: path to store acquired data to
    :param instrument_methods: methods and paths
    :param lib_id: defines the compound library
    :param method_suffix: defines the method, e.g., MSn, IT, HCD, ...
    :param plates_in_autosampler_location: list of tuples plate_id, location as tuples ("15", "B"),
    :param unique_id_header: defines a unique sample id. Must not end or start with a number so that contains matches are unique even for A1 and A10
    :param plate_id_header: plate id column. Plate ids can be any string or number
    :param well_header: defines the column with well locations, e.g., A1
    :param inject_volume_mul: micro liter injection volume
    :return:
    """
    # current_date = date.today().strftime("%Y%m%d")
    current_date = "20230620"
    # NO NEED TO CHANGE ANYTHING BELOW
    # final values
    data_filepath = os.path.join(data_filepath, lib_id, f"{current_date}_{method_suffix}")
    dataframes = []
    metadata_df = load_metadata_df(metadata_file, well_header, unique_id_header, method_suffix)
    for plate_id, plate_location in plates_in_autosampler_location:
        plate_df = filter_metadata_by_plate_id(metadata_df, plate_id, plate_id_header)

        sequence_file = f"data/Sequence/{current_date}_seq_rack_{plate_location}_{lib_id}_{plate_id}_{method_suffix}"

        df = _create_orbitrap_sequence(plate_df, sequence_file, data_filepath, well_header, plate_location,
                                       instrument_methods, inject_volume_mul)
        dataframes.append(df)
    concat = pd.concat(dataframes)
    plates_str = "_".join(["{}in{}".format(plate_id, loc) for plate_id, loc in plates_in_autosampler_location])
    sequence_file = f"data/Sequence/{current_date}_{plates_str}_seq_combined.csv"
    # add blanks and qcs
    final_df = add_blank_qc_rows(concat, data_filepath, instrument_methods, blank_well, qc_well,
                                 blank_qc_autosampler_location, blank_every_n_samples)

    write_thermo_sequence(sequence_file, final_df)


def add_blank_qc_rows(df: pd.DataFrame, data_filepath, instrument_methods, blank_well, qc_well,
                      blank_qc_autosampler_location, blank_every_n_samples) -> pd.DataFrame:
    if blank_well is None and qc_well is None:
        return df

    # TODO handle blank and qc autosampler location is None and use current plate

    main_method = instrument_methods[0]
    row_blank = {
        "File Name": "Blank",
        "Path": f"{data_filepath}_{main_method.identifier}",
        "Instrument Method": main_method.path,
        "Position": "{}:{}".format(blank_qc_autosampler_location, blank_well),
        "Inj Vol": "1",
        "Dil Factor": 1
    }
    row_qc = {
        "File Name": "QC",
        "Path": f"{data_filepath}_{main_method.identifier}",
        "Instrument Method": main_method.path,
        "Position": "{}:{}".format(blank_qc_autosampler_location, qc_well),
        "Inj Vol": "1",
        "Dil Factor": 1
    }
    blank_qc_df = pd.DataFrame([
        row_blank, row_qc
    ])
    chunks = pu.divide_chunks(df, blank_every_n_samples)
    chunks = [pd.concat([blank_qc_df, chunk]) for chunk in chunks]
    chunks.append(blank_qc_df)
    final_df = pd.concat(chunks)
    return final_df


def _create_orbitrap_sequence(metadata_df: DataFrame, sequence_file, data_filepath, well_header, plate_location,
                              instrument_methods: list[InstrumentMethod],
                              inject_volume_mul=3) -> DataFrame:
    """
    Creates Orbitrap sequence for positive and negative mode

    :param data_filepath: data acquisition file path. polarity will be added, path needs to be available before acquisition
            (date_suffix_polarity)
    :param sequence_file: the base sequence file to export. polarity and file type csv will be added automatically
    :param metadata_df: a dataframe that is already filtered to only contain samples from a single plate
    :param well_header: getting the well number of the final plate, e.g., A1
    :param instrument_methods: instrument methods
    :param plate_location: position in the autosampler
    Sequence for Plate 2 and Position Green
    :param inject_volume_mul: injection volume in micro liter
    :return:
    """
    if len(instrument_methods) == 0:
        raise ValueError("Provide at least one method file")

    dataframes = []
    for method in instrument_methods:
        polarity, instrument_method = method.identifier, method.path
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

    # current_date = date.today().strftime("%Y%m%d")
    current_date = "20230620"
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
