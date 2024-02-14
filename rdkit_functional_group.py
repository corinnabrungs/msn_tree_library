import rdkit.Chem as Chem
from dataclasses import dataclass, field
from enum import Enum
import pandas_utils as pu
import numpy as np
from smarts_utils import combine_smarts


@dataclass
class FunctionalGroup:
    group: str
    smarts: str
    pattern: Chem.rdchem.Mol = field(init=False)

    def __post_init__(self):
        self.pattern = Chem.MolFromSmarts(self.smarts)


def combine(groups):
    return combine_smarts([gr[1] for gr in groups])


class FunctionalGroups(FunctionalGroup, Enum):
    hydroxy = ("hydroxy", "[OX2H]")
    hydroxy_aliphatic = ("hydroxy_aliphatic", "[CX4][OX2H]")
    hydroxy_aromatic = ("hydroxy_aromatic", "[OX2H]c")
    sulfuric_acid_and_ester = (
        "sulfuric_acid_and_ester",
        "[$([SX4](=O)(=O)(O)O),$([SX4+2]([O-])([O-])(O)O)]",
    )
    sulfate = (
        "sulfate",
        "[$([#16X4](=[OX1])(=[OX1])([OX2H,OX1H0-])[OX2][#6]),$([#16X4+2]([OX1-])([OX1-])([OX2H,OX1H0-])[OX2][#6])]",
    )
    sulfuric_acid_diester = (
        "sulfuric_acid_diester",
        "[$([#16X4](=[OX1])(=[OX1])([OX2][#6])[OX2][#6]),$([#16X4](=[OX1])(=[OX1])([OX2][#6])[OX2][#6])]",
    )
    carbonyl = ("carbonyl", "[CX3]=[OX1]")
    acyl_halide = ("acyl_halide", "[CX3](=[OX1])[F,Cl,Br,I]")
    aldehyde = ("aldehyde", "[CX3H1](=O)[#6]")
    ketone = ("ketone", "[#6][CX3](=O)[#6]")
    carboxylic_acid = ("carboxylic_acid", "[CX3](=O)[OX2H,OX1H0-]")
    ester = ("ester", "[#6][CX3](=O)[OX2H0][#6]")
    lactone = ("lactone", "[$([O;R]@[C;R](=O))]")
    anhydride = ("anhydride", "[CX3](=[OX1])[OX2][CX3](=[OX1])")
    amide = ("amide", "[#7X3][#6X3](=O)")
    prim_amide = ("prim_amide", "[#7X3H2][#6X3](=[OX1])[#6]")
    second_amide = ("second_amide", "[#6][#7X3H1][#6X3](=[OX1])[#6]")
    tert_amide = ("tert_amide", "[#6][#7X3]([#6])[#6X3](=[OX1])[#6]")
    general_amide = ("general_amide", "[NX3][CX3](=[OX1])[#6]")
    amidinium = ("amidinium", "[NX3][CX3]=[NX3+]")
    carbamate = ("carbamate", "[NX3,NX4+][CX3](=[OX1])[OX2,OX1-]")
    cyanamide = ("cyanamide", "[NX3][CX2]#[NX1]")
    prim_amine = ("prim_amine", "[NX3;H2;!$([#7]C=[!#6]);!$([#7]C#[!#6])][#6]")
    second_amine = ("second_amine", "[#7X3H1;!$([#7][#6]=O)]")
    tert_amine = ("tert_amine", "[#7X3H0;!$([#7][#6]=O)]")
    quart_amine = ("quart_amine", "[#7X4H0+;!$([#7][#6]=O)]")
    lactam = ("lactam", "[$([N;R]@[C;R](=O))]")
    enamine = ("enamine", "[NX3][CX3]=[CX3]")
    aniline = ("aniline", "[NH2X3](cc)")
    prim_aromatic_amine = ("aromatic_amine", "[NX3H2][$(c)]")
    amino_acid = ("amino_acid", "[NX3,NX4+][CX4H]([*])[CX3](=[OX1])[O,N]")
    guanidine = (
        "guanidine",
        "[#6X3]([#7])(:,=[#7X2])[#7X3]",  # allow aromatic or double bond
    )
    azo_nitrogen = ("azo_nitrogen", "[NX2]=N")
    diazo_nitrogen = ("diazo_nitrogen", "[$([#6]=[N+]=[N-]),$([#6-]-[N+]#[N])]")
    azole = ("azole", "[$([nr5]:[nr5,or5,sr5]),$([nr5]:[cr5]:[nr5,or5,sr5])]")
    hydrazine = ("hydrazine", "[NX3][NX3]")
    hydrazone = ("hydrazone", "[NX3][NX2]=[*]")
    imine_schiff_base = (
        "imine_schiff_base",
        "[CX3;$([C]([#6])[#6]),$([CH][#6])]=[NX2][#6]",
    )
    nitrile = ("nitrile", "[NX1]#[CX2]")
    nitro = ("nitro", "[$([NX3](=O)=O),$([NX3+](=O)[O-])][!#8]")
    nitroso = ("nitroso", "[NX2]=[OX1]")
    enole = ("enole", "[OX2H][#6X3]=[#6]")
    phosphoric_acid = (
        "phosphoric_acid",
        "[$(P(=[OX1])([$([OX2H]),$([OX1-]),$([OX2]P)])([$([OX2H]),$([OX1-]),$([OX2]P)])[$([OX2H]),$([OX1-]),$([OX2]P)]),$([P+]([OX1-])([$([OX2H]),$([OX1-]),$([OX2]P)])([$([OX2H]),$([OX1-]),$([OX2]P)])[$([OX2H]),$([OX1-]),$([OX2]P)])]",
    )
    phosphoric_ester = (
        "phosphoric_ester",
        "[$(P(=[OX1])([OX2][#6])([$([OX2H]),$([OX1-]),$([OX2][#6])])[$([OX2H]),$([OX1-]),$([OX2][#6]),$([OX2]P)]),$([P+]([OX1-])([OX2][#6])([$([OX2H]),$([OX1-]),$([OX2][#6])])[$([OX2H]),$([OX1-]),$([OX2][#6]),$([OX2]P)])]",
    )
    carbo_thioester = ("carbo_thioester", "S([#6])[CX3](=O)[#6]")
    sulfoxide = ("sulfoxide", "[$([#16X3]=[OX1]),$([#16X3+][OX1-])]")
    sulfone = ("sulfone", "[$([#16X4](=[OX1])=[OX1]),$([#16X4+2]([OX1-])[OX1-])]")
    sulfonic_ester = (
        "sulfonic_ester",
        "[$([#16X4](=[OX1])(=[OX1])([#6])[OX2H0]),$([#16X4+2]([OX1-])([OX1-])([#6])[OX2H0])]",
    )
    sulfonic_acid = (
        "sulfonic_acid",
        "[$([#16X4](=[OX1])(=[OX1])([#6])[OX2H,OX1H0-]),$([#16X4+2]([OX1-])([OX1-])([#6])[OX2H,OX1H0-])]",
    )
    hexose = (
        "hexose",
        "[$(C(C1C(C(C(C([O]-1)[OX2])[OX2])[OX2])[OX2])[OX2]),$(C(C1C(C(C([O]-1)(C[OX2])[OX2])[OX2])[OX2])[OX2])]",
    )
    deoxy_hexose = (
        "deoxy_hexose",
        "[C;!$([#6][#8])](C1C(C(C(C([O]-1)[OX2])[OX2])[OX2])[OX2])",
    )
    pentose = (
        "pentose",
        "[$([C;!$([C][C])]1C(C(C(C([O]-1)[OX2])[OX2])[OX2])[OX2]),$(C(C1(C(C([C;!$([C][C])]([O]-1)[OX2])[OX2])[OX2]))[OX2])]",
    )
    glycoside = (
        "glycoside",
        combine([hexose, pentose, deoxy_hexose]),
    )
    flavan = (
        "flavan",
        "[#6]1-[#6]-[#6]2:[#6]:[#6]:[#6]:[#6]:[#6]:2-[#8]-[#6]-1-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1",
    )
    isoflavan = (
        "isoflavan",
        "[#6]1-[#6](-[#6]-[#8]-[#6]2:[#6]:[#6]:[#6]:[#6]:[#6]:2-1)-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1",
    )
    flavone = (
        "flavone",
        "[#6]1:[#6]:[#6]:[#6](:[#6]:[#6]:1)-[#6]1:[#6]:[#6](=[#8]):[#6]2:[#6]:[#6]:[#6]:[#6]:[#6]:2:[#8]:1",
    )
    isoflavone = (
        "isoflavone",
        "[#6]1:[#6]:[#6]:[#6](:[#6]:[#6]:1)-[#6]1:[#6]:[#8]:[#6]2:[#6]:[#6]:[#6]:[#6]:[#6]:2:[#6]:1=[#8]",
    )
    steroid = ("steroid", "[#6]1-[#6]-[#6]-[#6]2-[#6](-[#6]-1)-[#6]-[#6]-[#6]1-[#6]-2-[#6]-[#6]-[#6]2-[#6]-1-[#6]-[#6]-[#6]-2")
    # name = ("name", "C(C1C(C(C(C([O]-1)[OX2])[OX2])[OX2])[OX2])[!#8]")


def count_functional_groups(
    df, mol_col, groups: list[FunctionalGroup] = FunctionalGroups
):
    for group in groups:
        col = f"fg_n_{group.group}"
        df[col] = [count_functional_group(group, mol) for mol in mol_col]
        df = pu.astype_int(df, col)
    return df


def count_functional_group(group, mol) -> int:
    if pu.notnull(mol):
        return len(mol.GetSubstructMatches(group.pattern))
    else:
        return np.NAN
