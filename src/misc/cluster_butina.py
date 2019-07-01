#! /usr/bin/env python
#
#

import sys
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs
from rdkit.ML.Cluster import Butina


def ClusterFps(fps, cutoff=0.2):
    # first generate the distance matrix:
    dists = []
    nfps = len(fps)
    for i in range(1, nfps):
        sims = DataStructs.BulkTanimotoSimilarity(fps[i], fps[:i])
        dists.extend([1 - x for x in sims])

    # now cluster the data:
    cs = Butina.ClusterData(dists, nfps, cutoff, isDistData=True)
    return cs


def get_mols_from_file(fname):
    smiles_list = []
    logS_list = []
    with open(fname, 'r') as fin:
        for line in fin:
            pairs = line.rstrip('\n').split(',')
            smiles = pairs[0]
            logS = pairs[1]
            smiles_list.append(smiles)
            logS_list.append(logS)

    mols = [Chem.MolFromSmiles(x) for x in smiles_list]
    return mols, smiles_list, logS_list


if __name__ == "__main__":
    ms, smiles_list, logS_list = get_mols_from_file(sys.argv[1])
    fps = [AllChem.GetMorganFingerprintAsBitVect(x, 2, 1024) for x in ms]
    clusters = ClusterFps(fps, cutoff=0.4)

    with open(sys.argv[2], 'w') as fout:
        for cluster in clusters:
            exemplar_id = cluster[0]
            smiles = smiles_list[exemplar_id]
            logS = logS_list[exemplar_id]
            if float(logS) > 0.5 or float(logS) < -10.0:
                continue
            fout.write("{},{}\n".format(smiles, logS))
