import cirpy
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Descriptors import MolWt


def mol_wt(smiles):
    """Get molecular weight (in Daltons)"""
    return MolWt(Chem.MolFromSmiles(smiles))


def canonicalize_smiles(smiles, sanitize=True, iso=False, SLN=False):
    """Canonicalize given SMILES string
    The function is a wrapper around RDKIT function

    :argumnts:
      smiles -- (string) a compound in SMILES format
      sanitize -- (bool) sanitize the molecule
      iso -- (bool) include isomeric data in SMILES
      SLN -- (bool) is the molecule given in SLN format

    :return:
      canonicalized SMILES
    """
    if SLN:
        smiles_ = cirpy.resolve(smiles, "smiles")

        mol = Chem.MolToSmiles(
            Chem.MolFromSmiles(smiles_), canonical=True, isomericSmiles=iso
        )
    else:
        mol = Chem.MolToSmiles(
            Chem.MolFromSmiles(smiles), canonical=True, isomericSmiles=iso
        )

    mol = Chem.MolFromSmiles(mol)
    if sanitize:
        mol.UpdatePropertyCache(strict=False)
        mol = Chem.RemoveHs(
            mol, implicitOnly=False, updateExplicitCount=True, sanitize=True
)
        Chem.SanitizeMol(
            mol, Chem.rdmolops.SanitizeFlags.SANITIZE_ALL, catchErrors=False
        )
        AllChem.AssignStereochemistry(
            mol, cleanIt=True, force=True, flagPossibleStereoCenters=True
        )
        return Chem.MolToSmiles(mol, canonical=True, isomericSmiles=iso)
    else:
        return Chem.MolToSmiles(mol, canonical=True, isomericSmiles=iso)


if __name__ == "__main__":

    smiles_1 = "C1=CC=CN=C1"
    smiles_2 = "c1cccnc1"
    smiles_3 = "n1ccccc1"
    smiles_4 = "C1=CC=C(C(=C1)CC(=O)[O-])NC2=C(C=CC=C2Cl)Cl.[Na+]"
    smiles_5 = "CC(C)NCC(COC1=CC=C(C=C1)CC(=O)N)O"

    print(canonicalize_smiles(smiles_1))
    print(canonicalize_smiles(smiles_2))
    print(canonicalize_smiles(smiles_3))
    print(canonicalize_smiles(smiles_4))
    print(canonicalize_smiles(smiles_5))

    print(canonicalize_smiles("O=CHCH2CH2CH2CH2CH2CH2CH2CH3", SLN=True))
