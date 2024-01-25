from rdkit.Chem import MolStandardize, rdMolDescriptors, MolFromSmiles
import rdkit.Chem as Chem


def harmonize_smiles_rdkit(smiles, tautomer_limit=900):
    try:
        # take the largest covalently bound molecule
        smiles_largest = MolStandardize.fragment.LargestFragmentChooser(
            smiles
        ).prefer_organic
        mol = Chem.MolFromSmiles(smiles_largest)
        monomass = rdMolDescriptors.CalcExactMolWt(mol)

        # standardize tautomer
        if monomass < tautomer_limit:
            smiles_largest = MolStandardize.canonicalize_tautomer_smiles(smiles_largest)
            mol = Chem.MolFromSmiles(smiles_largest)

        # remove unnecessary charges
        uc = MolStandardize.charge.Uncharger()
        uncharged_mol = uc.uncharge(mol)

        # standardize the molecule
        lfc = MolStandardize.fragment.LargestFragmentChooser()
        standard_mol = lfc.choose(uncharged_mol)

        # remove stereochemistry
        Chem.RemoveStereochemistry(standard_mol)

        # get the standardized SMILES
        standard_smiles = Chem.MolToSmiles(standard_mol)
        return standard_smiles
    except Exception as e:
        print(f"An error occurred with input {smiles}: {e}")
        return ""
