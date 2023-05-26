from unittest import TestCase

import pandas as pd
from pubchem_client import pubchem_search_structure_by_name


class Test(TestCase):
    def test_pubchem_search_structure_by_name(self):
        df = pd.DataFrame(
            {
                "smiles": [None, "lfjlajdfi", "dawd", "dawdawd"],
                "compound_name": ["aspirin", "piperine", None, ""]
            }
        )
        df = pubchem_search_structure_by_name(df)

        assert len(df) == 5
