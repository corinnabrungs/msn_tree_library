import pandas as pd
from datetime import date
from tqdm import tqdm
import logging

tqdm.pandas()
logging.basicConfig(format="%(asctime)s - %(message)s", level=logging.DEBUG)


def main():
    metadata_file = r"data/nih/nih_library.tsv"

    # define all variables
    lib_id = "pluskal_nih"
    suffix = "_IT"  # is added to the end of the path and file names
    plate_id_header = "mixed_location_plate"
    current_date = date.today().strftime("%Y%m%d")

    instrument_method_positive = r"C:\Xcalibur\methods\Corinna_Brungs\IT_acquisition\IT_100AG_MS5_POS_mz115-2000"
    instrument_method_negative = r"C:\Xcalibur\methods\Corinna_Brungs\IT_acquisition\IT_100AG_MS5_NEG_mz115-2000"

    # plates are inserted into x compartment

    plates = ["15P"]
    plate_loc_in_autosampler = ["B"]
    # final values
    unique_id_header = "lib_plate_well"
    well_header = "final_plate_location"
    dataframes = []

    for filter_plate, plate_location in zip(plates, plate_loc_in_autosampler):
        df = create_orbitrap_sequence(lib_id, suffix, metadata_file, plate_id_header, unique_id_header, well_header,
                                      instrument_method_positive, instrument_method_negative, plate_location,
                                      filter_plate)
        dataframes.append(df)

    concat = pd.concat(dataframes)

    plates_str = "_".join(["{}in{}".format(plate, loc) for plate, loc in zip(plates, plate_loc_in_autosampler)])
    csv_file = f"data/Sequence/{current_date}_{plates_str}_seq_combined.csv"
    write_thermo_sequence(csv_file, concat)


def create_orbitrap_sequence(lib_id, suffix, metadata_file, plate_id_header, unique_id_header, well_header,
                             instrument_method_positive, instrument_method_negative, plate_location,
                             filter_plate) -> pd.DataFrame:
    """
    Creates Orbitrap sequence for positive and negative mode

    :param lib_id: name of your library/analysis
    :param suffix: suffix is added to the sequence file name right after the polarity - without an underscore _ so add one if wanted
    :param metadata_file: file path
    :param plate_id_header: file column header of plate and well location, e.g., Plate1_WellA1
    :param unique_id_header: combining the lib_id + plate + well
    :param well_header: getting the well number of the final plate, e.g., A1
    :param instrument_method_positive: instrument method in positive mode
    :param instrument_method_negative: instrument method in negative mode
    :param plate_location: position in the autosampler
    :param filter_plate: for multiple plates and position in the autosampler, Sequence for Plate 1 and Position Blue
    Sequence for Plate 2 and Position Green
    :return:
    """
    if suffix is None:
        suffix = ''

    current_date = date.today().strftime("%Y%m%d")
    logging.info("Will run on %s", metadata_file)

    instrument_methods = [instrument_method_positive, instrument_method_negative]
    polarities = ["positive", "negative"]
    dataframes = []
    for polarity, instrument_method in zip(polarities, instrument_methods):
        # import df
        if metadata_file.endswith(".tsv"):
            df = pd.read_csv(metadata_file, sep="\t")
        else:
            df = pd.read_csv(metadata_file, sep=",")

        # _id is added after the unique ID to securely substring search for the ID with number at the end
        # otherwise A1_1 is also contained A1_10 --> needs suffix
        df[unique_id_header] = ["{}_{}_id".format(lib_id, plate_id) for plate_id in df[plate_id_header]]
        df[well_header] = [plate_id.split("_")[1] for plate_id in df[plate_id_header]]
        seq_df = pd.DataFrame()
        seq_df["File Name"] = ["{}_{}_{}{}".format(current_date, unique_id, polarity, suffix) for unique_id in
                               df[unique_id_header]]
        #
        seq_df["Path"] = r"C:\Xcalibur\data\Corinna_Brungs\{}\{}_{}{}".format(lib_id, current_date, polarity, suffix)
        seq_df["Instrument Method"] = instrument_method
        seq_df["Position"] = ["{}:{}".format(plate_location, well) for well in df[well_header]]
        seq_df["Inj Vol"] = 3
        seq_df["Dil Factor"] = 1
        seq_df = seq_df.drop_duplicates()

        filtered_df = seq_df[seq_df["File Name"].str.contains(filter_plate)]
        csv_file = "data/Sequence/{}_seq_rack_{}_{}_{}_{}{}.csv".format(current_date, plate_location, lib_id,
                                                                        filter_plate, polarity, add)
        write_thermo_sequence(csv_file, filtered_df)
        dataframes.append(filtered_df)
    return pd.concat(dataframes)


def write_thermo_sequence(csv_file, df: pd.DataFrame):
    df.to_csv(csv_file, index=False)
    # df.to_csv("data/nih/{}_uniqueID.csv".format(lib_id), index=False)
    # Adding the first line as needed by Xcalibur sequence
    with open(csv_file, 'r') as original:
        data = original.read()
    with open(csv_file, 'w') as modified:
        modified.write("Bracket Type=4,\n" + data)


if __name__ == "__main__":
    main()
