import logging
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors, AllChem as Chem
from rdkit.Chem.MolStandardize import rdMolStandardize

from chembl_structure_pipeline import standardizer
from tqdm import tqdm

from meta_constants import MetaColumns
from pandas_utils import notnull, isnull, remove_empty_strings
import rdkit_functional_group
import rdkit_atom_count


# returns canonical smiles
def mol_to_canon_smiles(mol):
    try:
        return Chem.MolToSmiles(mol, isomericSmiles=False)
    except:
        return None


# returns isomerical smiles (if information is available)
def mol_to_isomeric_smiles(mol):
    try:
        return Chem.MolToSmiles(mol, isomericSmiles=True)
    except:
        return None


def mol_to_smarts(mol):
    try:
        return Chem.MolToSmarts(mol)
    except:
        return None


def chembl_standardize_mol(mol):
    return standardizer.standardize_mol(standardizer.get_parent_mol(mol)[0])


def smiles_to_mol(smiles: str):
    uncharger = rdMolStandardize.Uncharger()
    # smiles_stats = {'n_dots': Counter(), 'charge': Counter(), 'invalid_smiles': []}
    original_input = smiles
    try:
        smiles = split_smiles_major_mol(smiles)

        mol = Chem.MolFromSmiles(smiles)
        charge = Chem.GetFormalCharge(mol)
        if abs(charge) > 0:
            # smiles_stats['charge'][charge] += 1
            mol = uncharger.uncharge(mol)

        # if mol is None:
        #     return mol_from_pepseq(original_input)
        else:
            return mol
    except:
        return None


def split_smiles_major_mol(smiles):
    """

    :param smiles:
    :return: largest smiles sub part after splitting at . (salts etc)
    """
    # find the longest smiles that might be the main molecule
    # for smiles that contain the salt partner etc
    split_smiles = str(smiles).split(".")
    if len(split_smiles) > 1:
        return max(split_smiles, key=len)
    else:
        return split_smiles[0]


def split_inchikey(inchikey: str | None) -> str | None:
    return str(inchikey).split("-")[0] if notnull(inchikey) else None


def exact_mass_from_mol(mol):
    try:
        # canonical
        return round(Descriptors.ExactMolWt(mol), 6)
    except:
        return None


def inchi_from_mol(mol):
    try:
        return Chem.MolToInchi(mol)
    except:
        return None


def inchikey_from_mol(mol):
    try:
        return Chem.MolToInchiKey(mol)
    except:
        return None


def formula_from_mol(mol):
    try:
        return Chem.CalcMolFormula(mol)
    except:
        return None


def get_rdkit_mol(smiles, inchi):
    mol = None
    try:
        mol = Chem.MolFromSmiles(smiles)
    except:
        pass
    if mol is None:
        try:
            mol = Chem.MolFromInchi(inchi)
        except:
            pass
    return mol


def clean_structure_add_mol_id_columns(df, drop_mol=True) -> pd.DataFrame:
    """
    Will be performed twice to make sure structures are cleaned
    :param df: data frame with smiles and/or inchi columns
    :param drop_mol: drop mol column if True otherwise retain the column as "mol"
    :return: dataframe
    """
    logging.info("RDkit - predict properties")
    df = _add_molid_columns(df)
    df = _add_molid_columns(df)
    if drop_mol:
        df = df.drop(columns=["mol"], errors="ignore")
    return df


def _add_molid_columns(df) -> pd.DataFrame:
    if MetaColumns.inchi not in df.columns:
        df[MetaColumns.inchi] = None
    # merge all smiles from isomeric_smiles>canonical_smiles>smiles
    df = ensure_smiles_column(df)
    # first strip any salts
    df[MetaColumns.smiles] = [
        split_smiles_major_mol(smiles) if notnull(smiles) else np.NAN
        for smiles in df["smiles"]
    ]
    df["mol"] = [
        get_rdkit_mol(smiles, inchi) for smiles, inchi in zip(df["smiles"], df["inchi"])
    ]
    df["mol"] = [
        chembl_standardize_mol(mol) if notnull(mol) else np.NAN
        for mol in tqdm(df["mol"], "clean structure")
    ]
    df[MetaColumns.canonical_smiles] = [mol_to_canon_smiles(mol) for mol in df["mol"]]
    df[MetaColumns.isomeric_smiles] = [mol_to_isomeric_smiles(mol) for mol in df["mol"]]
    df[MetaColumns.smarts] = [mol_to_smarts(mol) for mol in df["mol"]]
    df[MetaColumns.monoisotopic_mass] = [exact_mass_from_mol(mol) for mol in df["mol"]]
    df[MetaColumns.inchi] = [inchi_from_mol(mol) for mol in df["mol"]]
    df[MetaColumns.inchikey] = [inchikey_from_mol(mol) for mol in df["mol"]]
    df[MetaColumns.split_inchikey] = [
        split_inchikey(inchikey) for inchikey in df["inchikey"]
    ]
    df[MetaColumns.formula] = [formula_from_mol(mol) for mol in df["mol"]]
    df[MetaColumns.logp] = [
        Descriptors.MolLogP(mol) if notnull(mol) else np.NAN for mol in df["mol"]
    ]
    # counting specific function groups
    rdkit_functional_group.count_functional_groups(df, df["mol"])
    rdkit_atom_count.count_element_atoms_df(df, df["mol"])

    # merge all smiles from isomeric_smiles>canonical_smiles>smiles
    df = ensure_smiles_column(df)
    return df


def ensure_smiles_column(df: pd.DataFrame) -> pd.DataFrame:
    """
    Ensure presence of smiles column, remove empty strings, and fill all none with smiles in
    isomeric_smiles>canonical_smiles>smiles
    """
    if "inchi" not in df.columns:
        df["inchi"] = None
    if "isomerical_smiles" in df.columns:
        df = df.rename(columns={"isomerical_smiles": MetaColumns.isomeric_smiles})

    # ensure smiles column by priority
    headers = [
        MetaColumns.isomeric_smiles,
        MetaColumns.canonical_smiles,
        MetaColumns.smiles,
        "SMILES",
        "Smiles",
    ]

    headers = [h for h in headers if h in df.columns]

    if len(headers)==0:
        df[MetaColumns.smiles] = None
        return df

    df = remove_empty_strings(df, headers)

    # keep columns temporarily
    cols = [df[h] for h in headers]
    df[MetaColumns.smiles] = cols[0]

    if len(headers)>1:
        # fill NA values by priority
        for col in cols[1:]:
            df[MetaColumns.smiles] = df[MetaColumns.smiles].fillna(col)
    return df
