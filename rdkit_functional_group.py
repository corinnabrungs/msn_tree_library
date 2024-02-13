import rdkit.Chem as Chem
from dataclasses import dataclass, field
import pandas_utils as pu
import numpy as np


@dataclass
class FunctionalGroup:
    name: str
    smarts: str
    pattern: Chem.rdchem.Mol = field(init=False)

    def __post_init__(self):
        self.pattern = Chem.MolFromSmarts(self.smarts)


hydroxy = FunctionalGroup("hydroxy", "[*OH]")
sulfuric_acid = FunctionalGroup(
    "sulfuric_acid", "[$([SX4](=O)(=O)(O)O),$([SX4+2]([O-])([O-])(O)O)]"
)
sulfate = FunctionalGroup(
    "sulfate",
    "[$([#16X4](=[OX1])(=[OX1])([OX2H,OX1H0-])[OX2][#6]),$([#16X4+2]([OX1-])([OX1-])([OX2H,OX1H0-])[OX2][#6])]",
)
# name = FunctionalGroup("hydroxy", "place")


default_groups = [hydroxy, sulfate, sulfuric_acid]


def count_functional_groups(
    df, mol_col, groups: list[FunctionalGroup] = default_groups
):
    for group in groups:
        col = f"fg_n_{group.name}"
        df[col] = [count_functional_group(group, mol) for mol in mol_col]
        df = pu.astype_int(df, col)
    return df


def count_functional_group(group, mol) -> int:
    if pu.notnull(mol):
        return len(mol.GetSubstructMatches(group.pattern))
    else:
        return np.NAN
