import pandas as pd
import logging
from drug_utils import map_clinical_phase_to_number
from tqdm import tqdm
from rdkit_mol_identifiers import split_inchikey

tqdm.pandas()


def broad_list_search(df):
    logging.info("Search broad institute list of drugs by first block of inchikey")
    # download from: https://clue.io/repurposing#download-data
    prefix = "broad_"
    broad_df = pd.read_csv("data/broad_institute_drug_list.csv")
    broad_df = broad_df.add_prefix(prefix)

    if "split_inchikey" not in df and "inchikey" in df:
        df["split_inchikey"] = [split_inchikey(inchikey) for inchikey in df['inchikey']]
    broad_df["split_inchikey"] = [split_inchikey(inchikey) for inchikey in broad_df["{}InChIKey".format(prefix)]]

    merged_df = pd.merge(df, broad_df, on="split_inchikey", how="left")
    # converting the clinical phases (from broad institute, chembl, provider, or else) to numbers (remove phase,
    # preclinic (as 0.5), or launched)
    merged_df["broad_clinical_phase"] = [map_clinical_phase_to_number(phase) for phase in
                                         merged_df["broad_clinical_phase"]]
    merged_df["clinical_phase"] = [map_clinical_phase_to_number(phase) for phase in merged_df["clinical_phase"]]
    if "Clinical Information" in df.columns:
        merged_df["Clinical Information"] = [map_clinical_phase_to_number(phase) for phase in
                                             merged_df["Clinical Information"]]
    else:
        merged_df["Clinical Information"] = None

    # Comparing the clinical phases, only store the highest number in clinical phase column
    merged_df["clinical_phase"] = merged_df[['broad_clinical_phase', 'clinical_phase', 'Clinical Information']].max(
        axis=1)
    return merged_df.drop(["broad_clinical_phase", "Clinical Information", "{}InChIKey".format(prefix)], axis=1)
