from unittest import TestCase

import drugcentral_client
import drugcentral_postgresql_query as drugcentral
from meta_constants import MetaColumns


class TestDrugCentral(TestCase):

    def test_search_split_inchikey(self):
        try:
            drugcentral.connect()
            columns, structures = drugcentral.drugcentral_postgresql(split_inchikey="LVWZTYCIRDMTEY")
        finally:
            drugcentral.deconnect()
        assert structures
        assert structures[0] == 4659
        assert structures[2] == 'metamizole'

    def test_database_not_connected(self):
        columns, structures = drugcentral.drugcentral_postgresql(inchikey="SWMDAPWAQQTBOG-UHFFFAOYSA-N")
        assert structures == None

    # def test_to_dataframe(self):
    #     import pandas as pd
    #     inchikeys = ["SWMDAPWAQQTBOG-UHFFFAOYSA-N", "SWMDAPWAQQTBOG-UHFFFAOYSA-N", "SWMDAPWAQQTBOG-UHFFFAOYSA-N"]
    #
    #     try:
    #         drugcentral.connect()
    #         results = [drugcentral.drugcentral_postgresql(inchikey=inchikey) for inchikey in inchikeys]
    #         data = [row[1] for row in results]
    #         first_entry_columns = next((row[0] for row in results), [])
    #         columns = [col.name for col in first_entry_columns]
    #         assert len(pd.DataFrame(data=data, columns=columns))==3
    #     finally:
    #         drugcentral.deconnect()

    def test_merge_df(self):
        import pandas as pd
        df = pd.DataFrame({
            "inchikey": ["GKDRMWXFWHEQQT-UHFFFAOYSA-N", "SWMDAPWAQQTBOG-WRONG-N", None],
            "madeup_name": ["Correct", "Wrong", "Nothing"],
        })
        merged_df = drugcentral_client.drugcentral_search(df)
        assert len(merged_df) == 3
        assert merged_df.at[0, "drugcentral_cas"] == '901119-35-5'
        assert merged_df.at[0, "drugcentral_id"] == '5280'
        assert merged_df.at[1, "drugcentral_compound_name"] == 'fostemsavir'
        assert merged_df[MetaColumns.date_drugcentral_search].isnull().sum() == 0

    def test_by_ids(self):
        import pandas as pd
        df = pd.DataFrame.from_dict({
            "identifier": ["", "inchikey", "chembl_id", "pubchem_cid", "unii", "drugbank_id", "split_inchikey"],
            "inchikey": ["", "GKDRMWXFWHEQQT-UHFFFAOYSA-N", "", "", "", "", ""],
            "chembl_id": ["", "", "CHEMBL298734", "", "", "", ""],
            "pubchem_cid": ["", "", "chembl", "68165256", "", "", ""],
            "unii": ["", "", "chembl", "", "P76B05O5V6", "", ""],
            "drugbank_id": ["", "", "chembl", "", "", "DB09252", ""],
            "split_inchikey": ["", "", "CHEMBL", "", "", "", "GKDRMWXFWHEQQT"],
        })
        merged_df = drugcentral_client.drugcentral_search(df)
        assert len(merged_df) == 7
        assert len(merged_df.columns) > 10
