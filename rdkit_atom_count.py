import pandas as pd
import rdkit.Chem as Chem
import pandas_utils as pu
import numpy as np

from dataclasses import dataclass, field
from enum import Enum


@dataclass
class Element:
    symbol: str
    atomic_number: int


class Elements(Element, Enum):
    H = ("H", 1)
    C = ("C", 6)
    N = ("N", 7)
    O = ("O", 8)
    P = ("P", 15)
    S = ("S", 16)
    F = ("F", 9)
    Cl = ("Cl", 17)
    Br = ("Br", 35)
    I = ("I", 53)


def count_element_atoms_df(
    df, mol_col, elements: list[Element] = Elements
) -> pd.DataFrame:
    # also count H
    mol_col_h = [Chem.AddHs(mol) for mol in mol_col]

    for element in elements:
        col = f"at_n_{element.symbol}"
        df[col] = [count_element(mol, element.atomic_number) for mol in mol_col_h]
        df = pu.astype_int(df, col)
    return df


def count_element(mol: Chem.rdchem.Mol, number: int) -> int:
    return sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == number)
