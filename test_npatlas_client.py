from unittest import TestCase
import pandas as pd

import npatlas_client


class Test(TestCase):
    def test_search_np_atlas(self):
        df = pd.DataFrame([
            {
                "compound_name": "piperine",
                "smiles": "C1CCN(CC1)C(=O)/C=C/C=C/C2=CC3=C(C=C2)OCO3",
                "inchikey": "MXXWOMGUGJBKIW-YPCIICBESA-N",
                "inchi": "InChI=1S/C17H19NO3/c19-17(18-10-4-1-5-11-18)7-3-2-6-14-8-9-15-16(12-14)21-13-20-15/h2-3,6-9,12H,1,4-5,10-11,13H2/b6-2+,7-3+",
            },
            {
                "compound_name": "piperine",
                "inchikey": "",
            },
            {
                "compound_name": "piperine",
                "inchikey": None,
            },
        ])
        df = npatlas_client.search_np_atlas(df)

        new_cols = [col for col in df.columns if npatlas_client.NP_ATLAS_PREFIX in col]
        assert len(new_cols) > 2
