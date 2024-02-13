from unittest import TestCase
import rdkit_functional_group as fg
import rdkit.Chem as Chem


def mol(smiles):
    return Chem.MolFromSmiles(smiles)


def count(group: fg.FunctionalGroup, smiles: str):
    return fg.count_functional_group(group, mol(smiles))


class Test(TestCase):
    def test_count(self):
        cases = [
            (fg.guanidine, "CS(=O)(=O)N1CCC2=C(N=C(N=C21)N3CCOCC3)C4=CN=C(N=C4)N", 2),
        ]

        for case in cases:
            group, smiles, expected = case
            found = count(group, smiles)
            self.assertEqual(
                found,
                expected,
                f"{group.name} smarts {group.smarts} failed count in {smiles}",
            )


if __name__ == "__main__":
    unittest.main()
