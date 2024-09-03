import pandas as pd
import ast
import pyteomics.mgf
import numpy as np
from tqdm.notebook import tqdm
from pandas_utils import isnull, notnull


def extract_energies(energy, msn_energies):
    if notnull(energy):
        return energy
    if isnull(msn_energies):
        return None

    return msn_energies[-1]


def read_mgf(infile, lib=None) -> pd.DataFrame:
    import re

    rows = []
    counter = 0

    with pyteomics.mgf.MGF(str(infile)) as f_in:
        for spectrum_dict in tqdm(f_in):
            if spectrum_dict is not None:
                rows.append(spectrum_dict["params"])
            else:
                counter += 1

    df = pd.DataFrame(rows)
    if "inchikey" not in df.columns:
        df["inchikey"] = df["inchiaux"]
    if "compound_name" not in df.columns:
        df["compound_name"] = df["name"]
    if "monoisotopic_mass" not in df.columns:
        df["monoisotopic_mass"] = df["exactmass"]
    if "spectype" not in df.columns:
        df["spectype"] = None
    df["spectype"] = df["spectype"].fillna("BEST_TIC")
    if "msn_collision_energies" in df.columns:
        df["msn_collision_energies"] = [
            ast.literal_eval(energies) if notnull(energies) else np.NAN
            for energies in (df["msn_collision_energies"])
        ]
        df["collision energy"] = [
            extract_energies(energy, msn_energies)
            for energy, msn_energies in zip(
                df["collision energy"], df["msn_collision_energies"]
            )
        ]
    if "quality_explained_intensity" in df.columns:
        df["quality_explained_intensity"] = df["quality_explained_intensity"].astype(
            float
        )
    if "quality_explained_signals" in df.columns:
        df["quality_explained_signals"] = df["quality_explained_signals"].astype(float)

    #
    if lib != None and "usi" in df.columns:
        df["unique_sample_id"] = [
            "{}{}_id".format(lib, re.search(rf"{lib}(.*?)_id", usi).group(1))
            for usi in df["usi"]
        ]
    return df
