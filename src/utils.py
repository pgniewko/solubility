
import cirpy
from rdkit import Chem
from rdkit.Chem.Descriptors import MolWt


def mol_wt(smiles):
    return MolWt(Chem.MolFromSmiles(smiles))


def canonicalize_smiles(smiles, iso=False, SLN=False):
    if SLN:
        smiles_ = cirpy.resolve(smiles, 'smiles')
        return Chem.MolToSmiles(Chem.MolFromSmiles(smiles_), canonical=True, isomericSmiles=iso)
    else:
        return Chem.MolToSmiles(Chem.MolFromSmiles(smiles), canonical=True, isomericSmiles=iso)


if __name__ == "__main__":

    smiles_1 = 'C1=CC=CN=C1'
    smiles_2 = 'c1cccnc1'
    smiles_3 = 'n1ccccc1'
    smiles_4 = 'C1=CC=C(C(=C1)CC(=O)[O-])NC2=C(C=CC=C2Cl)Cl.[Na+]'
    smiles_5 = 'CC(C)NCC(COC1=CC=C(C=C1)CC(=O)N)O'

    print(canonicalize_smiles(smiles_1))
    print(canonicalize_smiles(smiles_2))
    print(canonicalize_smiles(smiles_3))
    print(canonicalize_smiles(smiles_4))
    print(canonicalize_smiles(smiles_5))

    print(canonicalize_smiles("O=CHCH2CH2CH2CH2CH2CH2CH2CH3", SLN=True))
