from unittest import TestCase
import pandas as pd

import lotus_client


class Test(TestCase):
    def test_apply_search_on_split_inchikey(self):
        df = pd.DataFrame(
            {
                "split_inchikey": ["XQYZDYMELSJDRZ", "XQYZDYMELSJDRZ", "BSYNRYMUTXBXSQ", "noinchikey", None, None]
            }
        )
        results = lotus_client.apply_search_on_split_inchikey(df)

        self.assertEqual(len(results.columns), 14)
        self.assertEqual(len(results), 6)
        assert results.at[0, "parent_taxon_rank_label"] == "genus; order"
        assert results.at[2, "taxon_name"] == "Glycyrrhiza glabra; Ixora coccinea"

    def test_read_lotus_dataframe(self):
        df = lotus_client.read_lotus_dataframe()
        self.assertIsNotNone(df, "Cannot import lotus data, run the prepare_wikidata_lotus_data_prefect.py script")
        self.assertEqual(len(df.columns), 14)
        self.assertEqual(len(df), 148738)
