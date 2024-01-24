from unittest import TestCase
import my_mces
import pandas as pd


class Test(TestCase):
    def test_mces_parallel(self):
        df = pd.DataFrame(
            {
                "smiles1": [
                    "NC(=O)c1c(NC(=O)C2CCCCC2)sc(Nc2ccccc2)n1",
                    "O=Cc1c[nH]c2c1c(O)ccc2",
                    "C[C@H](O)[C@H](O)c1cnc2[nH]c(N)nc(=O)c2n1",
                    "OC[C@H]1O[C@@](CO)(OC[C@@]2(O[C@H]3O[C@H](CO)[C@@H](O)[C@H](O)[C@H]3O)O[C@H](CO)[C@@H](O)[C@@H]2O)[C@@H](O)[C@@H]1O",
                ],
                "smiles2": [
                    "OCC(O)C(O)C(OC1OC(CO)C(O)C(O)C1O)C(O)CO",
                    "Oc1nc2c(O)cccc2cc1",
                    "CC(O)C(O)c1cnc2[nH]c(=N)nc(O)c2n1",
                    "OC[C@H]1O[C@H](OC[C@H]2O[C@H](O[C@]3(CO)O[C@H](CO)[C@@H](O)[C@@H]3O)[C@H](O)[C@@H](O)[C@@H]2O)[C@H](O)[C@@H](O)[C@H]1O",
                ],
            }
        )

        my_mces.mces_parallel(df, "smiles1", "smiles2", 4, "mces_results.tsv")

    def test_ceil(self):
        from math import ceil

        print(ceil(3 / 2))
