from unittest import TestCase
import pandas as pd

import structure_classifier_client
from rest_utils import get_json_response
from structure_classifier_client import np_classifier_url


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

    def test_np_classifier(self):
        smiles = "COc1c(O)ccc(/C=C/C(=O)NCC(O)c2ccc(O)cc2)c1"
        npclass = np_classifier_url(smiles)
        json = get_json_response(npclass)

        assert npclass != None
        assert json != None
