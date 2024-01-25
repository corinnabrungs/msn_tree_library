from unittest import TestCase
import structure_cleanup_rdkit as clean


class Test(TestCase):
  def test_harmonize_smiles_rdkit(self):
    smiles = clean.harmonize_smiles_rdkit("O=S(C1=C2C=C(S(=O)(O[Na])=O)C(/N=N/C3=C4C=C(S(=O)(O[Na])=O)C=CC4=C(/N=N/C5=CC=C(S(=O)(O[Na])=O)C=C5)C=C3)=C(O)C2=C(NC(C)=O)C=C1)(O[Na])=O")
    print(smiles)