import logging

import numpy as np
import pandas as pd
from chembl_webresource_client.new_client import new_client as chembl

from date_utils import create_expired_entries_dataframe, iso_datetime_now
from meta_constants import MetaColumns
from pandas_utils import (
    notnull,
    isnull,
    update_dataframes,
    make_str_floor_to_int_number,
)
from drug_utils import map_clinical_phase_to_number
from tqdm import tqdm
import datetime as dt


def get_chembl_mol(chembl_id=None, inchikey=None):
    try:
        if isnull(chembl_id) and isnull(inchikey):
            return None

        if notnull(chembl_id):
            comp = chembl.molecule.get(chembl_id)
            if comp:
                return comp

        compounds = None
        if not compounds and notnull(inchikey):
            compounds = chembl.molecule.filter(
                molecule_structures__standard_inchi_key=inchikey
            )
        if not compounds:
            logging.info(
                "NO ChEMBL FOR: chemblid: {} or inchikey: {}".format(
                    chembl_id, inchikey
                )
            )
            return None
        else:
            return compounds[0]
    except Exception as e:
        logging.warning("Error during chembl query:", e)
        return None


def chembl_search_id_and_inchikey(
    df, refresh_expired_entries_after: dt.timedelta = dt.timedelta(days=90)
) -> pd.DataFrame:
    logging.info("Search ChEMBL by chemblid or inchikey")
    if "chembl_id" not in df.columns:
        df["chembl_id"] = None

    # only work on expired elements
    # define which rows are old or were not searched before
    filtered = create_expired_entries_dataframe(
        df, MetaColumns.date_chembl_search, refresh_expired_entries_after
    )

    filtered["result_column"] = [
        get_chembl_mol(chembl_id, inchikey)
        for chembl_id, inchikey in tqdm(
            zip(filtered["chembl_id"], filtered["inchikey"])
        )
    ]
    filtered = filtered[filtered["result_column"].notnull()].copy()
    compounds = filtered["result_column"]
    # refresh date
    filtered[MetaColumns.date_chembl_search] = iso_datetime_now()

    filtered["chembl_id"] = [compound["molecule_chembl_id"] for compound in compounds]
    # filtered["compound_name"] = filtered["compound_name"] + [compound["pref_name"] for compound in compounds]
    filtered["prodrug"] = [compound["prodrug"] for compound in compounds]
    filtered["availability"] = [compound["availability_type"] for compound in compounds]
    filtered["chembl_clinical_phase"] = [
        compound["max_phase"] for compound in compounds
    ]
    filtered["chembl_clinical_phase"] = [
        map_clinical_phase_to_number(phase)
        for phase in filtered["chembl_clinical_phase"]
    ]
    filtered["withdrawn"] = [compound["withdrawn_flag"] for compound in compounds]
    filtered[MetaColumns.first_approval] = pd.array(
        [compound["first_approval"] for compound in compounds], dtype=pd.Int64Dtype()
    )
    filtered = make_str_floor_to_int_number(filtered, MetaColumns.first_approval)
    filtered["oral"] = [compound["oral"] for compound in compounds]
    filtered["parenteral"] = [compound["parenteral"] for compound in compounds]
    filtered["topical"] = [compound["topical"] for compound in compounds]
    filtered["natural_product"] = [
        compound["natural_product"] for compound in compounds
    ]
    filtered["usan_stem_definition"] = [
        compound["usan_stem_definition"] for compound in compounds
    ]
    filtered["chembl_indication"] = [
        compound["indication_class"] for compound in compounds
    ]
    filtered["chembl_atc_classifications"] = [
        compound["atc_classifications"] for compound in compounds
    ]

    # properties sometimes None
    props = [compound["molecule_properties"] for compound in compounds]

    filtered["molecular_species"] = [
        prop["molecular_species"] if notnull(prop) else None for prop in props
    ]
    filtered["chembl_alogp"] = [
        prop["alogp"] if notnull(prop) else None for prop in props
    ]
    filtered["chembl_cx_logp"] = [
        prop["cx_logp"] if notnull(prop) else None for prop in props
    ]

    # was changed by ChEMBL api
    # filtered["withdrawn_class"] = [compound["withdrawn_class"] for compound in compounds]
    # filtered["withdrawn_reason"] = [compound["withdrawn_reason"] for compound in compounds]
    # filtered["withdrawn_year"] = pd.array([compound["withdrawn_year"] for compound in compounds], dtype=pd.Int64Dtype())
    # filtered["withdrawn_country"] = [compound["withdrawn_country"] for compound in compounds]

    ## dont overwrite
    # filtered["synonyms"] = filtered["synonyms"] + [compound["molecule_synonyms"] if notnull(compound) else [] for
    # compound in compounds]

    # combine new data with old rows that were not processed
    return update_dataframes(filtered, df).drop(
        columns=["result_column"], errors="ignore"
    )
