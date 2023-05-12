import logging

import numpy as np
import pandas as pd
from chembl_webresource_client.new_client import new_client as chembl
from tqdm import tqdm


def get_chembl_mol(chembl_id=None, inchikey=None):
    try:
        if not (chembl_id or inchikey):
            raise ValueError("At least one chembl identifier needs to be a value")

        if chembl_id:
            comp = chembl.molecule.get(chembl_id)
            if comp:
                return comp

        compounds = None
        if not compounds and inchikey:
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

    df["chembl_id"] = [compound["molecule_chembl_id"] if pd.notnull(compound) else np.NAN for compound in compounds]
    # df["compound_name"] = df["compound_name"] + [compound["pref_name"] if pd.notnull(compound) else np.NAN for
    # compound in compounds]
    df["molecular_species"] = [
        compound["molecule_properties"]["molecular_species"] if pd.notnull(compound) else np.NAN for compound in
        compounds]
    df["prodrug"] = [compound["prodrug"] if pd.notnull(compound) else np.NAN for compound in compounds]
    df["availability"] = [compound["availability_type"] if pd.notnull(compound) else np.NAN for compound in compounds]
    df["clinical_phase"] = [compound["max_phase"] if pd.notnull(compound) else np.NAN for compound in compounds]
    df["first_approval"] = pd.array(
        [compound["first_approval"] if pd.notnull(compound) else np.NAN for compound in compounds],
        dtype=pd.Int64Dtype())
    df["withdrawn"] = [compound["withdrawn_flag"] if pd.notnull(compound) else np.NAN for compound in compounds]
    # was changed by ChEMBL api
    # df["withdrawn_class"] = [compound["withdrawn_class"] if pd.notnull(compound) else np.NAN for compound in
    #                          compounds]
    # df["withdrawn_reason"] = [compound["withdrawn_reason"] if pd.notnull(compound) else np.NAN for compound in
    #                           compounds]
    # df["withdrawn_year"] = pd.array(
    #     [compound["withdrawn_year"] if pd.notnull(compound) else np.NAN for compound in compounds],
    #     dtype=pd.Int64Dtype())
    # df["withdrawn_country"] = [compound["withdrawn_country"] if pd.notnull(compound) else np.NAN for compound in
    #                            compounds]
    df["oral"] = [compound["oral"] if pd.notnull(compound) else np.NAN for compound in compounds]
    df["parenteral"] = [compound["parenteral"] if pd.notnull(compound) else np.NAN for compound in compounds]
    df["topical"] = [compound["topical"] if pd.notnull(compound) else np.NAN for compound in compounds]
    df["natural_product"] = [compound["natural_product"] if pd.notnull(compound) else np.NAN for compound in
                             compounds]
    df["usan_stem_definition"] = [compound["usan_stem_definition"] if pd.notnull(compound) else np.NAN for compound
                                  in compounds]
    df["chembl_alogp"] = [compound["molecule_properties"]["alogp"] if pd.notnull(compound) else np.NAN for compound
                          in compounds]
    df["chembl_clogp"] = [compound["molecule_properties"]["cx_logp"] if pd.notnull(compound) else np.NAN for compound
                          in compounds]

    ## dont overwrite
    # df["synonyms"] = df["synonyms"] + [compound["molecule_synonyms"] if pd.notnull(compound) else [] for
    # compound in compounds]
    df["indication"] = [compound["indication_class"] if pd.notnull(compound) else np.NAN for compound in compounds]

    return df
