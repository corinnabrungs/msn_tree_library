from unittest import TestCase
import pandas as pd

import broadinstitute_client


class Test(TestCase):
    def test_get_broadinstitute_columns(self):
        df = pd.DataFrame(
            [
                {
                    "inchikey": "XDLYKKIQACFMJG-UHFFFAOYSA-N",
                    "smiles": None,
                    "clinical_phase": "1",
                }
            ]
        )
        df = broadinstitute_client.broad_list_search(df)
        print(df)
