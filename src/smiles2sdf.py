#! /usr/bin/env python

import sys
import numpy as np

from rdkit import Chem


def read_smiles(fname):
    smiles_list = []
    with open(fname, 'r') as fin:
        for line in fin:
            smiles = line.rstrip('\n').split(',')[0]
            smiles_list.append(smiles)
   

    return [Chem.MolFromSmiles(x) for x in smiles_list]


if __name__ == "__main__":
    fin = sys.argv[1]
    fout = sys.argv[2]
   
    mols = read_smiles(fin)
    w = Chem.SDWriter(fout)
    for m in mols: 
        w.write(m)
