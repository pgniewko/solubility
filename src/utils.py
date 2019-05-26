from rdkit import Chem
from standardiser import standardise

def canonicalize_smiles(smiles, iso=False):
    return Chem.MolToSmiles(Chem.MolFromSmiles(smiles), canonical=True, isomericSmiles=iso)

def prepare_smiles(smiles_orig, iso=False):
    mol = Chem.MolFromSmiles(smiles_orig)
    msg = None
    try:
        Chem.SanitizeMol(mol, Chem.rdmolops.SanitizeFlags.SANITIZE_ALL, False)
        noh_mol = Chem.RemoveHs(mol)
    except ValueError as e:
        msg = "Failed sanitizning molecule: %s" %(str(e))
        return None, msg

    try:
        clean_mol = standardise.run(noh_mol)
    except standardise.StandardiseException as e:
        msg = "Failed to standardize molecule %s" %(str(e))
        return None, msg
    except RuntimeError:
        msg = "Failed to parse SMILES"
        return None, msg

    smiles = Chem.MolToSmiles(clean_mol, canonical=True, isomericSmiles=iso)
    return smiles, None



if __name__ == "__main__":
    smiles_1 = 'C1=CC=CN=C1'
    smiles_2 = 'c1cccnc1'
    smiles_3 = 'n1ccccc1'
    
    print (canonicalize_smiles(smiles_1))
    print (canonicalize_smiles(smiles_2))
    print (canonicalize_smiles(smiles_3))
     
    print (prepare_smiles(smiles_1))
    print (prepare_smiles(smiles_2))
    print (prepare_smiles(smiles_3))
