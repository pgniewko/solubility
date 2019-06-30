#! /usr/bin/env python

import sys
from rdkit import Chem


def read_smiles(fname):
    smiles_list = []
    logS = []
    with open(fname, 'r') as fin:
        for line in fin:
            pairs = line.rstrip('\n').split(',')
            smiles = pairs[0]
            if len(pairs) > 1:
                logS.append(float(pairs[1]))

            smiles_list.append(smiles)

    return [Chem.MolFromSmiles(x) for x in smiles_list], logS


if __name__ == "__main__":
    fin = sys.argv[1]
    fout = sys.argv[2]

    mols, logS = read_smiles(fin)
    w = Chem.SDWriter(fout)
    for i, m in enumerate(mols):
        try:
            m.SetProp('logS', str(logS[i]))
        except IndexError:
            pass

        w.write(m)
