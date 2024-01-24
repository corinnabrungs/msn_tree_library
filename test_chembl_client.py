from unittest import TestCase
import pandas as pd

import chembl_client


class Test(TestCase):
    def test_get_chembl_mol(self):
        comp = chembl_client.get_chembl_mol(
            None, inchikey="BSYNRYMUTXBXSQ-UHFFFAOYSA-N"
        )

        assert comp is not None
        assert comp["max_phase"] == "4.0"
        assert comp["oral"] == True
        assert comp["pref_name"] == "ASPIRIN"

    def test_chembl_search_id_and_inchikey(self):
        df = pd.DataFrame(
            [
                {
                    "smiles": "CC(=O)Oc1ccccc1C(=O)O",
                    "inchikey": "BSYNRYMUTXBXSQ-UHFFFAOYSA-N",
                    "inchi": "InChI=1S/C9H8O4/c1-6(10)13-8-5-3-2-4-7(8)9(11)12/h2-5H,1H3,(H,11,12)",
                }
            ]
        )
        df = chembl_client.chembl_search_id_and_inchikey(df)

        cols = ["chembl_alogp", "indication", "oral", "first_approval", "prodrug"]
        for col in cols:
            assert col in df.columns
