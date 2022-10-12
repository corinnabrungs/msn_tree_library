from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import AllChem as Chem
from rdkit.Chem.MolStandardize import rdMolStandardize

from chembl_structure_pipeline import standardizer

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
    split_smiles = str(smiles).split('.')
    if len(split_smiles) > 1:
        return max(split_smiles, key=len)
    else:
        return split_smiles[0]


def exact_mass_from_mol(mol):
    try:
        # canonical
        return Descriptors.ExactMolWt(mol)
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