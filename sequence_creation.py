import pandas as pd
from datetime import date
from tqdm import tqdm
import logging

tqdm.pandas()
logging.basicConfig(format="%(asctime)s - %(message)s", level=logging.DEBUG)


def main():
    metadata_file = r"C:\git\msn_library\data\mce_library.tsv"

    # define all variables
    lib_id = "pluskal_mce"
    plate_id_header = "mixed_location_plate1"

    instrument_method_positive = r"C:\Xcalibur\methods\Corinna_Brungs\Library6_100AGC_60000Res_MS5_POS_mz115-2000"
    instrument_method_negative = r"C:\Xcalibur\methods\Corinna_Brungs\Library6_100AGC_60000Res_MS5_NEG_mz115-2000"

    # plates are inserted into the BLUE B compartment
    plates = ["1D1", "1D2", "1D3"]
    plate_loc_in_autosampler = ["B", "G", "R"]

    # final values
    unique_id_header = "lib_plate_well"
    well_header = "final_plate_location"

    for filter_plate, plate_location in zip(plates, plate_loc_in_autosampler):
        create_orbitrap_sequence(lib_id, metadata_file, plate_id_header, unique_id_header, well_header,
                                 instrument_method_positive, instrument_method_negative, plate_location, filter_plate)


def create_orbitrap_sequence(lib_id, metadata_file, plate_id_header, unique_id_header, well_header,
                             instrument_method_positive, instrument_method_negative, plate_location, filter_plate):
    """
    Creates Orbitrap sequence for positive and negative mode

    :param lib_id: name of your library/analysis
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
    current_date = date.today().strftime("%Y%m%d")
    logging.info("Will run on %s", metadata_file)

    instrument_methods = [instrument_method_positive, instrument_method_negative]
    polarities = ["positive", "negative"]
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
        seq_df["File Name"] = ["{}_{}_{}".format(current_date, unique_id, polarity) for unique_id in
                               df[unique_id_header]]
        seq_df["Path"] = r"C:\Xcalibur\data\Corinna_Brungs\{}".format(lib_id)
        seq_df["Instrument Method"] = instrument_method
        seq_df["Position"] = ["{}:{}".format(plate_location, well) for well in df[well_header]]
        seq_df["Inj Vol"] = 2
        seq_df["Dil Factor"] = 1
        seq_df = seq_df.drop_duplicates()

        filtered_df = seq_df[seq_df["File Name"].str.contains(filter_plate)]
        csv_file = "data/seq_rack_{}_{}_{}_{}.csv".format(plate_location, lib_id, filter_plate, polarity)
        filtered_df.to_csv(csv_file, index=False)

        # Adding the first line as needed by Xcalibur sequence
        with open(csv_file, 'r') as original:
            data = original.read()
        with open(csv_file, 'w') as modified:
            modified.write("Bracket Type=4,\n" + data)


if __name__ == "__main__":
    main()
