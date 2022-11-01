from unittest import TestCase

import drugcentral_postgresql_query as drugcentral
import metadata_cleanup


class TestDrugCentral(TestCase):

    def test_search_split_inchikey(self):
        try:
            drugcentral.connect()
            columns, structures = drugcentral.drugcentral_postgresql(split_inchikey="LVWZTYCIRDMTEY")
        finally:
            drugcentral.deconnect()
        assert structures


    def test_database_not_connected(self):
        columns, structures = drugcentral.drugcentral_postgresql(inchikey="SWMDAPWAQQTBOG-UHFFFAOYSA-N")
        assert structures == None


    def test_to_dataframe(self):
        import pandas as pd
        inchikeys = ["SWMDAPWAQQTBOG-UHFFFAOYSA-N", "SWMDAPWAQQTBOG-UHFFFAOYSA-N", "SWMDAPWAQQTBOG-UHFFFAOYSA-N"]

        try:
            drugcentral.connect()
            results = [drugcentral.drugcentral_postgresql(inchikey=inchikey) for inchikey in inchikeys]
            data = [row[1] for row in results]
            first_entry_columns = next((row[0] for row in results), [])
            columns = [col.name for col in first_entry_columns]
            assert len(pd.DataFrame(data=data, columns=columns))==3
        finally:
            drugcentral.deconnect()


    def test_merge_df(self):
        import pandas as pd
        df = pd.DataFrame({
            "inchi_key" : ["GKDRMWXFWHEQQT-UHFFFAOYSA-N", "SWMDAPWAQQTBOG-WRONG-N", None],
            "madeup_name" : ["Correct", "Wrong", "Nothing"],
        })
        merged_df = metadata_cleanup.drugcentral_search(df)
        assert len(merged_df) == 3