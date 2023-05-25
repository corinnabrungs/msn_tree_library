import logging

import numpy as np
import pandas as pd
from chembl_webresource_client.new_client import new_client as chembl
from pandas_utils import notnull, isnull
from tqdm import tqdm


def get_chembl_mol(chembl_id=None, inchikey=None):
    try:
        if isnull(chembl_id) and isnull(inchikey):
            raise ValueError("At least one chembl identifier needs to be a value")

        if notnull(chembl_id):
            comp = chembl.molecule.get(chembl_id)
            if comp:
                return comp

        compounds = None
        if not compounds and notnull(inchikey):
            compounds = chembl.molecule.filter(molecule_structures__standard_inchi_key=inchikey)
        if not compounds:
            logging.info("NO ChEMBL FOR: chemblid: {} or inchikey: {}".format(chembl_id, inchikey))
            return None
        else:
            return compounds[0]
    except Exception as e:
        logging.warning("Error during chembl query:", e)
        return None


def chembl_search_id_and_inchikey(df) -> pd.DataFrame:
    logging.info("Search ChEMBL by chemblid or inchikey")
    compounds = [get_chembl_mol(chembl_id, inchikey) for chembl_id, inchikey in
                 tqdm(zip(df["chembl_id"], df["inchikey"]))]

    df["chembl_id"] = [compound["molecule_chembl_id"] if notnull(compound) else np.NAN for compound in compounds]
    # df["compound_name"] = df["compound_name"] + [compound["pref_name"] if notnull(compound) else np.NAN for
    # compound in compounds]
    df["molecular_species"] = [
        compound["molecule_properties"]["molecular_species"] if notnull(compound) else np.NAN for compound in
        compounds]
    df["prodrug"] = [compound["prodrug"] if notnull(compound) else np.NAN for compound in compounds]
    df["availability"] = [compound["availability_type"] if notnull(compound) else np.NAN for compound in compounds]
    df["clinical_phase"] = [compound["max_phase"] if notnull(compound) else np.NAN for compound in compounds]
    df["first_approval"] = pd.array(
        [compound["first_approval"] if notnull(compound) else np.NAN for compound in compounds],
        dtype=pd.Int64Dtype())
    df["withdrawn"] = [compound["withdrawn_flag"] if notnull(compound) else np.NAN for compound in compounds]
    # was changed by ChEMBL api
    # df["withdrawn_class"] = [compound["withdrawn_class"] if notnull(compound) else np.NAN for compound in
    #                          compounds]
    # df["withdrawn_reason"] = [compound["withdrawn_reason"] if notnull(compound) else np.NAN for compound in
    #                           compounds]
    # df["withdrawn_year"] = pd.array(
    #     [compound["withdrawn_year"] if notnull(compound) else np.NAN for compound in compounds],
    #     dtype=pd.Int64Dtype())
    # df["withdrawn_country"] = [compound["withdrawn_country"] if notnull(compound) else np.NAN for compound in
    #                            compounds]
    df["oral"] = [compound["oral"] if notnull(compound) else np.NAN for compound in compounds]
    df["parenteral"] = [compound["parenteral"] if notnull(compound) else np.NAN for compound in compounds]
    df["topical"] = [compound["topical"] if notnull(compound) else np.NAN for compound in compounds]
    df["natural_product"] = [compound["natural_product"] if notnull(compound) else np.NAN for compound in
                             compounds]
    df["usan_stem_definition"] = [compound["usan_stem_definition"] if notnull(compound) else np.NAN for compound
                                  in compounds]
    df["chembl_alogp"] = [compound["molecule_properties"]["alogp"] if notnull(compound) else np.NAN for compound
                          in compounds]
    df["chembl_clogp"] = [compound["molecule_properties"]["cx_logp"] if notnull(compound) else np.NAN for compound
                          in compounds]

    ## dont overwrite
    # df["synonyms"] = df["synonyms"] + [compound["molecule_synonyms"] if notnull(compound) else [] for
    # compound in compounds]
    df["indication"] = [compound["indication_class"] if notnull(compound) else np.NAN for compound in compounds]

    return df
