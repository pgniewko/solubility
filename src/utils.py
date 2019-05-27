
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem.Descriptors import MolWt, ExactMolWt


def mol_wt(smiles):
    return MolWt(Chem.MolFromSmiles(smiles))

def canonicalize_smiles(smiles, iso=False):
    return Chem.MolToSmiles(Chem.MolFromSmiles(smiles), canonical=True, isomericSmiles=iso)

if __name__ == "__main__":
    smiles_1 = 'C1=CC=CN=C1'
    smiles_2 = 'c1cccnc1'
    smiles_3 = 'n1ccccc1'
    smiles_4 = 'C1=CC=C(C(=C1)CC(=O)[O-])NC2=C(C=CC=C2Cl)Cl.[Na+]'
    smiles_5 = 'CC(C)NCC(COC1=CC=C(C=C1)CC(=O)N)O'
    smiles_6 = 'O=C(CO)[C@@]3(O)CC[C@H]2[C@@H]4CC\C1=C\C(=O)CC[C@]1(C)[C@H]4C(=O)C[C@@]23C'
    
    print (canonicalize_smiles(smiles_1))
    print (canonicalize_smiles(smiles_2))
    print (canonicalize_smiles(smiles_3))
    print(canonicalize_smiles(smiles_4))
    print(canonicalize_smiles(smiles_5))
    print(canonicalize_smiles(smiles_6))
