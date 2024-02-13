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


hydroxy = FunctionalGroup("hydroxy", "[OX2H]")
hydroxy_aliphatic = FunctionalGroup("hydroxy_aliphatic", "[CX4][OX2H]")
phenole = FunctionalGroup("phenole", "[OX2H]c")
sulfuric_acid_and_ester = FunctionalGroup(
    "sulfuric_acid_and_ester", "[$([SX4](=O)(=O)(O)O),$([SX4+2]([O-])([O-])(O)O)]"
)
sulfate = FunctionalGroup(
    "sulfate",
    "[$([#16X4](=[OX1])(=[OX1])([OX2H,OX1H0-])[OX2][#6]),$([#16X4+2]([OX1-])([OX1-])([OX2H,OX1H0-])[OX2][#6])]",
)
sulfuric_acid_diester = FunctionalGroup(
    "sulfuric_acid_diester",
    "[$([#16X4](=[OX1])(=[OX1])([OX2][#6])[OX2][#6]),$([#16X4](=[OX1])(=[OX1])([OX2][#6])[OX2][#6])]",
)
carbonyl = FunctionalGroup("carbonyl", "[CX3]=[OX1]")
acyl_halide = FunctionalGroup("acyl_halide", "[CX3](=[OX1])[F,Cl,Br,I]")
aldehyde = FunctionalGroup("aldehyde", "[CX3H1](=O)[#6]")
ketone = FunctionalGroup("ketone", "[#6][CX3](=O)[#6]")
carboxylic_acid = FunctionalGroup("carboxylic_acid", "[CX3](=O)[OX2H,OX1H0-]")
ester = FunctionalGroup("ester", "[#6][CX3](=O)[OX2H0][#6]")
lactone = FunctionalGroup("lactone", "[$([O;R]@[C;R](=O))]")
anhydride = FunctionalGroup("anhydride", "[CX3](=[OX1])[OX2][CX3](=[OX1])")
amide = FunctionalGroup("amide", "[#7X3][#6X3](=O)")
prim_amide = FunctionalGroup("prim_amide", "[#7X3H2][#6X3](=[OX1])[#6]")
second_amide = FunctionalGroup("second_amide", "[#6][#7X3H1][#6X3](=[OX1])[#6]")
tert_amide = FunctionalGroup("tert_amide", "[#6][#7X3]([#6])[#6X3](=[OX1])[#6]")
general_amide = FunctionalGroup("general_amide", "[NX3][CX3](=[OX1])[#6]")
amidinium = FunctionalGroup("amidinium", "[NX3][CX3]=[NX3+]")
carbamate = FunctionalGroup("carbamate", "[NX3,NX4+][CX3](=[OX1])[OX2,OX1-]")
cyanamide = FunctionalGroup("cyanamide", "[NX3][CX2]#[NX1]")
prim_amine = FunctionalGroup(
    "prim_amine", "[NX3;H2;!$([#7]C=[!#6]);!$([#7]C#[!#6])][#6]"
)
second_amine = FunctionalGroup("second_amine", "[#7X3H1;!$([#7][#6]=O)]")
tert_amine = FunctionalGroup("tert_amine", "[#7X3H0;!$([#7][#6]=O)]")
quart_amine = FunctionalGroup("quart_amine", "[#7X4H0+;!$([#7][#6]=O)]")
enamine = FunctionalGroup("enamine", "[NX3][CX3]=[CX3]")
aniline = FunctionalGroup("aniline", "[NH2X3](cc)")
prim_aromatic_amine = FunctionalGroup("aromatic_amine", "[NX3H2][$(c)]")
amino_acid = FunctionalGroup("amino_acid", "[NX3,NX4+][CX4H]([*])[CX3](=[OX1])[O,N]")
guanidine = FunctionalGroup(
    "guanidine", "[#6X3]([#7])(:,=[#7X2])[#7X3]"  # allow aromatic or double bond
)
azo_nitrogen = FunctionalGroup("azo_nitrogen", "[NX2]=N")
diazo_nitrogen = FunctionalGroup(
    "diazo_nitrogen", "[$([#6]=[N+]=[N-]),$([#6-]-[N+]#[N])]"
)
azole = FunctionalGroup(
    "azole", "[$([nr5]:[nr5,or5,sr5]),$([nr5]:[cr5]:[nr5,or5,sr5])]"
)
hydrazine = FunctionalGroup("hydrazine", "[NX3][NX3]")
hydrazone = FunctionalGroup("hydrazone", "[NX3][NX2]=[*]")
imine_schiff_base = FunctionalGroup(
    "imine_schiff_base", "[CX3;$([C]([#6])[#6]),$([CH][#6])]=[NX2][#6]"
)
nitrile = FunctionalGroup("nitrile", "[NX1]#[CX2]")
nitro = FunctionalGroup("nitro", "[$([NX3](=O)=O),$([NX3+](=O)[O-])][!#8]")
nitroso = FunctionalGroup("nitroso", "[NX2]=[OX1]")
enole = FunctionalGroup("enole", "[OX2H][#6X3]=[#6]")
phosphoric_acid = FunctionalGroup(
    "phosphoric_acid",
    "[$(P(=[OX1])([$([OX2H]),$([OX1-]),$([OX2]P)])([$([OX2H]),$([OX1-]),$([OX2]P)])[$([OX2H]),$([OX1-]),$([OX2]P)]),$([P+]([OX1-])([$([OX2H]),$([OX1-]),$([OX2]P)])([$([OX2H]),$([OX1-]),$([OX2]P)])[$([OX2H]),$([OX1-]),$([OX2]P)])]",
)
phosphoric_ester = FunctionalGroup(
    "phosphoric_ester",
    "[$(P(=[OX1])([OX2][#6])([$([OX2H]),$([OX1-]),$([OX2][#6])])[$([OX2H]),$([OX1-]),$([OX2][#6]),$([OX2]P)]),$([P+]([OX1-])([OX2][#6])([$([OX2H]),$([OX1-]),$([OX2][#6])])[$([OX2H]),$([OX1-]),$([OX2][#6]),$([OX2]P)])]",
)
carbo_thioester = FunctionalGroup("carbo_thioester", "S([#6])[CX3](=O)[#6]")
sulfoxide = FunctionalGroup("sulfoxide", "[$([#16X3]=[OX1]),$([#16X3+][OX1-])]")
sulfone = FunctionalGroup(
    "sulfone", "[$([#16X4](=[OX1])=[OX1]),$([#16X4+2]([OX1-])[OX1-])]"
)
sulfonic_ester = FunctionalGroup(
    "sulfonic_ester",
    "[$([#16X4](=[OX1])(=[OX1])([#6])[OX2H0]),$([#16X4+2]([OX1-])([OX1-])([#6])[OX2H0])]",
)
sulfonic_acid = FunctionalGroup(
    "sulfonic_acid",
    "[$([#16X4](=[OX1])(=[OX1])([#6])[OX2H,OX1H0-]),$([#16X4+2]([OX1-])([OX1-])([#6])[OX2H,OX1H0-])]",
)
# name = FunctionalGroup("name", "place")
# name = FunctionalGroup("name", "place")
# name = FunctionalGroup("name", "place")


default_groups = [
    hydroxy,
    hydroxy_aliphatic,
    sulfate,
    sulfuric_acid_and_ester,
    sulfuric_acid_diester,
    carbonyl,
    acyl_halide,
    aldehyde,
    anhydride,
    amide,
    prim_amide,
    second_amide,
    tert_amide,
    general_amide,
    amidinium,
    carbamate,
    carboxylic_acid,
    cyanamide,
    ester,
    ketone,
    prim_amine,
    second_amine,
    tert_amine,
    quart_amine,
    enamine,
    aniline,
    prim_aromatic_amine,
    amino_acid,
    guanidine,
    azo_nitrogen,
    diazo_nitrogen,
    azole,
    hydrazine,
    hydrazone,
    imine_schiff_base,
    nitrile,
    nitro,
    nitroso,
    enole,
    phosphoric_ester,
    phosphoric_acid,
    carbo_thioester,
    sulfoxide,
    sulfone,
    sulfonic_ester,
    sulfonic_acid,
]


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
