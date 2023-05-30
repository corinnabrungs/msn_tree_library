from unittest import TestCase
import pandas as pd

import structure_classifier_client


class Test(TestCase):

    def test_apply_classyfire(self):
        df = pd.DataFrame([{
            "smiles": "CC(=O)Oc1ccccc1C(=O)O",
            "inchikey": "BSYNRYMUTXBXSQ-UHFFFAOYSA-N",
            "inchi": "InChI=1S/C9H8O4/c1-6(10)13-8-5-3-2-4-7(8)9(11)12/h2-5H,1H3,(H,11,12)",
        }])
        df = structure_classifier_client.apply_classyfire(df)

        new_cols = [col for col in df.columns if structure_classifier_client.CLASSYFIRE_PREFIX in col]
        assert len(new_cols) > 2

    def test_apply_np_classifier(self):
        df = pd.DataFrame([{
            "smiles": "CC(=O)Oc1ccccc1C(=O)O",
            "inchikey": "BSYNRYMUTXBXSQ-UHFFFAOYSA-N",
            "inchi": "InChI=1S/C9H8O4/c1-6(10)13-8-5-3-2-4-7(8)9(11)12/h2-5H,1H3,(H,11,12)",
        }])
        df = structure_classifier_client.apply_np_classifier(df)

        new_cols = [col for col in df.columns if structure_classifier_client.NP_CLASSIFIER_PREFIX in col]
        assert len(new_cols) > 2
