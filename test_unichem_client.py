from unittest import TestCase

import unichem_client
import pandas as pd


class Test(TestCase):
    def test_search_unichem_xref(self):
        inchikey = "MXXWOMGUGJBKIW-YPCIICBESA-N"
        result = unichem_client.search_unichem_xref(inchikey, search_type="inchikey")
        assert result is not None

    def test_search_unichem_xref_listids(self):
        inchikey = "XQYZDYMELSJDRZ-UHFFFAOYSA-N"
        result = unichem_client.search_unichem_xref(inchikey, search_type="inchikey")
        assert result is not None

    def test_extract_ids_to_columns(self):
        df = pd.DataFrame(
            {
                "inchikey": [
                    "XQYZDYMELSJDRZ-UHFFFAOYSA-N",
                    "XQYZDYMELSJDRZ-UHFFFAOYSA-N",
                    "XQYZDYMELSJDRZ-UHFFFAOYSA-N",
                    "MXXWOMGUGJBKIW-YPCIICBESA-N",
                    "BSYNRYMUTXBXSQ-UHFFFAOYSA-N",
                    None,
                    "",
                ]
            }
        )
        df = unichem_client.search_all_xrefs(df)

        assert unichem_client.Sources.unichem_source.column_name in df.columns
        assert unichem_client.Sources.pubchem_source.column_name in df.columns
        assert unichem_client.Sources.chembl_source.column_name in df.columns
