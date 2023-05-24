from unittest import TestCase
import pandas as pd

import drugbank_client


class Test(TestCase):
    def test_drugbank_search_add_columns(self):
        df = pd.DataFrame(
            [
                {
                    "inchikey": "MXXWOMGUGJBKIW-YPCIICBESA-N",
                    "drugbank_id": "DB12582"
                },
                {
                    "inchikey": None,
                    "drugbank_id": "DB00945"
                },
                {
                    "inchikey": "MXXWOMGUGJBKIW-YPCIICBESA-N",
                    "drugbank_id": None
                }
            ]
        )

        result = drugbank_client.drugbank_search_add_columns(df)

        assert len(result) > 1
        assert result.at[0, "compound_name"] == "Piperine"
        assert result.at[1, "chembl_id"] == "CHEMBL25"
        assert result.at[1, "chembl_id"] != "CHEMBL26"
