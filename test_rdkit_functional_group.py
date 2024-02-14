import unittest
from unittest import TestCase
import rdkit.Chem as Chem
from dataclasses import dataclass

from rdkit_functional_group import count_functional_group, FunctionalGroup
from rdkit_functional_group import FunctionalGroups as fg


def mol(smiles):
    return Chem.MolFromSmiles(smiles)


def count(group: FunctionalGroup, smiles: str):
    return count_functional_group(group, mol(smiles))


@dataclass
class Case:
    expected_matches: int
    smiles: str
    group: FunctionalGroup


class Test(TestCase):
    def test_count(self):
        taurocholicacid = (
            "CC(CCC(=O)NCCS(=O)(=O)O)C1CCC2C1(C(CC3C2C(CC4C3(CCC(C4)O)C)O)O)C"
        )
        CH5132799 = "CS(=O)(=O)N1CCC2=C(N=C(N=C21)N3CCOCC3)C4=CN=C(N=C4)N"
        fostemsavir_trom = "CC1=NN(C=N1)C2=NC=C(C3=C2N(C=C3C(=O)C(=O)N4CCN(CC4)C(=O)C5=CC=CC=C5)COP(=O)(O)O)OC.C(C(CO)(CO)N)O"
        phenol = "C1=CC=C(C=C1)O"
        ascorbicacid = "C(C(C1C(=C(C(=O)O1)O)O)O)O"
        patulin = "C1C=C2C(=CC(=O)O2)C(O1)O"
        folic_acid = "C1=CC(=CC=C1C(=O)NC(CCC(=O)O)C(=O)O)NCC2=CN=C3C(=N2)C(=O)NC(=N3)N"
        folate = "C1=CC(=CC=C1C(=O)NC(CCC(=O)[O-])C(=O)O)NCC2=CN=C3C(=N2)C(=O)NC(=N3)N"
        lecitinpc = "CCCCCCCCCCCCCCCC(=O)OCC(COP(=O)([O-])OCC[N+](C)(C)C)OC(=O)CCCCCCCC=CCC=CCCCCC"
        penicillin_g = "CC1(C(N2C(S1)C(C2=O)NC(=O)CC3=CC=CC=C3)C(=O)O)C"
        poliumoside = "CC1C(C(C(C(O1)OCC2C(C(C(C(O2)OCCC3=CC(=C(C=C3)O)O)O)OC4C(C(C(C(O4)C)O)O)O)OC(=O)C=CC5=CC(=C(C=C5)O)O)O)O)O"
        rutin = "CC1C(C(C(C(O1)OCC2C(C(C(C(O2)OC3=C(OC4=CC(=CC(=C4C3=O)O)O)C5=CC(=C(C=C5)O)O)O)O)O)O)O)O"

        cases = [
            Case(2, CH5132799, fg.guanidine),
            Case(0, CH5132799, fg.sulfonic_ester),
            Case(0, CH5132799, fg.sulfonic_acid),
            Case(1, CH5132799, fg.sulfone),
            Case(0, CH5132799, fg.sulfoxide),
            Case(0, CH5132799, fg.sulfate),
            Case(1, CH5132799, fg.prim_aromatic_amine),
            Case(1, phenol, fg.hydroxy_aromatic),
            # taurocholic acid
            Case(0, taurocholicacid, fg.hydroxy_aromatic),
            Case(4, taurocholicacid, fg.hydroxy),
            Case(3, taurocholicacid, fg.hydroxy_aliphatic),
            Case(1, taurocholicacid, fg.amide),
            Case(0, taurocholicacid, fg.sulfonic_ester),
            Case(1, taurocholicacid, fg.sulfonic_acid),
            Case(0, taurocholicacid, fg.lactone),
            Case(0, fostemsavir_trom, fg.sulfonic_acid),
            Case(0, fostemsavir_trom, fg.prim_aromatic_amine),
            Case(1, fostemsavir_trom, fg.prim_amine),
            Case(0, fostemsavir_trom, fg.second_amine),
            Case(2, fostemsavir_trom, fg.tert_amine),
            Case(0, fostemsavir_trom, fg.quart_amine),
            Case(0, fostemsavir_trom, fg.enamine),
            Case(0, fostemsavir_trom, fg.enole),
            Case(5, fostemsavir_trom, fg.hydroxy),
            Case(3, fostemsavir_trom, fg.hydroxy_aliphatic),
            Case(1, fostemsavir_trom, fg.phosphoric_acid),
            Case(1, fostemsavir_trom, fg.phosphoric_ester),
            Case(0, fostemsavir_trom, fg.guanidine),
            Case(2, fostemsavir_trom, fg.amide),
            Case(2, fostemsavir_trom, fg.tert_amide),
            Case(0, fostemsavir_trom, fg.prim_amide),
            Case(0, fostemsavir_trom, fg.second_amide),
            Case(3, fostemsavir_trom, fg.azole),
            Case(1, fostemsavir_trom, fg.ketone),
            Case(0, fostemsavir_trom, fg.aldehyde),
            Case(3, fostemsavir_trom, fg.carbonyl),
            Case(0, fostemsavir_trom, fg.carboxylic_acid),
            Case(0, fostemsavir_trom, fg.lactone),
            Case(0, ascorbicacid, fg.carboxylic_acid),
            Case(0, ascorbicacid, fg.ketone),
            Case(1, ascorbicacid, fg.carbonyl),
            Case(1, ascorbicacid, fg.lactone),
            Case(4, ascorbicacid, fg.hydroxy),
            Case(1, patulin, fg.hydroxy),
            Case(1, patulin, fg.lactone),
            Case(1, patulin, fg.ester),
            Case(0, folic_acid, fg.lactone),
            Case(0, folic_acid, fg.ester),
            Case(2, folic_acid, fg.amide),
            Case(2, folic_acid, fg.second_amide),
            Case(1, folic_acid, fg.prim_amine),
            Case(1, folic_acid, fg.second_amine),
            Case(0, folic_acid, fg.tert_amine),
            Case(2, folic_acid, fg.carboxylic_acid),
            Case(0, folate, fg.lactone),
            Case(0, folate, fg.ester),
            Case(2, folate, fg.amide),
            Case(2, folate, fg.second_amide),
            Case(1, folate, fg.prim_amine),
            Case(1, folate, fg.second_amine),
            Case(0, folate, fg.tert_amine),
            Case(2, folate, fg.carboxylic_acid),
            Case(0, lecitinpc, fg.carboxylic_acid),
            Case(0, lecitinpc, fg.prim_amine),
            Case(0, lecitinpc, fg.second_amine),
            Case(0, lecitinpc, fg.tert_amine),
            Case(1, lecitinpc, fg.quart_amine),
            Case(1, lecitinpc, fg.phosphoric_acid),
            Case(1, lecitinpc, fg.phosphoric_ester),
            Case(0, lecitinpc, fg.lactone),
            Case(2, lecitinpc, fg.ester),
            Case(0, penicillin_g, fg.lactone),
            Case(1, penicillin_g, fg.lactam),
            Case(1, penicillin_g, fg.carboxylic_acid),
            Case(2, penicillin_g, fg.amide),
            Case(0, penicillin_g, fg.prim_amide),
            Case(1, penicillin_g, fg.second_amide),
            Case(1, penicillin_g, fg.tert_amide),
            Case(2, penicillin_g, fg.amino_acid),
            Case(0, penicillin_g, fg.ketone),
            Case(0, poliumoside, fg.ketone),
            Case(11, poliumoside, fg.hydroxy),
            Case(7, poliumoside, fg.hydroxy_aliphatic),
            Case(4, poliumoside, fg.hydroxy_aromatic),
            Case(1, poliumoside, fg.hexose),
            Case(2, poliumoside, fg.deoxy_hexose),
            Case(3, poliumoside, fg.glycoside),
            Case(0, poliumoside, fg.pentose),
            Case(0, "C(C1C(C(C(O1)(CO)O)O)O)O", fg.pentose),
            Case(1, "C(C1C(C(C(O1)(CO)O)O)O)O", fg.hexose),
            Case(0, "C(C1C(C(C(O1)(CO)O)O)O)O", fg.deoxy_hexose),
            Case(0, "C(C1C(C(C(O1)(CO)O)O)O)O", fg.deoxy_hexose),
            Case(1, rutin, fg.deoxy_hexose),
            Case(1, rutin, fg.hexose),
            Case(2, rutin, fg.glycoside),
            Case(0, rutin, fg.flavan),
            Case(1, rutin, fg.flavone),
        ]

        for case in cases:
            group = case.group
            found = count(group, case.smiles)
            self.assertEqual(
                case.expected_matches,
                found,
                f"{group.group} smarts {group.smarts} failed count in {case.smiles}",
            )


if __name__ == "__main__":
    unittest.main()
