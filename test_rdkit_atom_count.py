import unittest
from unittest import TestCase

import pandas as pd

import rdkit_atom_count as ac
import rdkit.Chem as Chem
from rdkit_atom_count import Elements


def count_atom(smiles, element: Elements):
    return ac.count_element(mol(smiles), element.atomic_number)


def mol(smiles):
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)
    return mol


class Test(TestCase):
    def test_count_element(self):
        phenol = "C1=CC=C(C=C1)O"
        self.check(6, Elements.C, phenol)
        self.check(1, Elements.O, phenol)
        self.check(6, Elements.H, phenol)

    def test_count_element_df(self):
        df = pd.DataFrame({"mol": [mol("C1=CC=C(C=C1)O")]})
        df = ac.count_element_atoms_df(df, df["mol"])
        self.assertEqual(6, df.at[0, f"at_n_{Elements.C.symbol}"])
        self.assertEqual(1, df.at[0, f"at_n_{Elements.O.symbol}"])
        self.assertEqual(6, df.at[0, f"at_n_{Elements.H.symbol}"])

    def test_count_element_df_limit_elements(self):
        df = pd.DataFrame({"mol": [mol("C1=CC=C(C=C1)O")]})
        # only search C
        df = ac.count_element_atoms_df(df, df["mol"], [Elements.C])
        self.assertEqual(6, df.at[0, f"at_n_{Elements.C.symbol}"])
        self.assertTrue(f"at_n_{Elements.C.symbol}" in df.columns)
        self.assertFalse(f"at_n_{Elements.H.symbol}" in df.columns)
        self.assertFalse(f"at_n_{Elements.O.symbol}" in df.columns)

    def check(self, expected, element: ac.Elements, smiles):
        self.assertEqual(expected, count_atom(smiles, element))


if __name__ == "__main__":
    unittest.main()
