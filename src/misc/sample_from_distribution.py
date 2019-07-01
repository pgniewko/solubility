#! /usr/bin/env python
# Generate uniform distribution from another distribution.
# 

import sys
import random
import numpy as np
from collections import defaultdict


def get_bins_probs(logS_list, min_logS, max_logS, dx):
    nbins = int((max_logS - min_logS) / dx)
    hist, edges = np.histogram(logS_list, bins=nbins, range=(min_logS, max_logS), normed=True)
    prob_norm = np.sum(1. / (dx * hist))
    probs = (1. / (dx * hist)) / prob_norm

    bins_dict = defaultdict(list)
    for i, val in enumerate(logS_list):
        if val < edges[0] or val > edges[-1]:
            continue
        val_shifted = val - edges[0]
        idx = int(val_shifted / dx)
        bins_dict[idx].append(i)

    return hist, edges, probs, bins_dict


def get_cmpd_idx(bins_dict):
    N = len(bins_dict.keys())
    bin_no = random.randint(0, N-2)
    vals = bins_dict[bin_no]
    return random.choice(vals)


if __name__ == "__main__":

    N = 10000
    logS_list = []
    smiles_list = []
    with open(sys.argv[1], 'r') as fin:
        for line in fin:
            pairs = line.rstrip('\n').split(',')
            smiles = pairs[0]
            logS = float(pairs[1])

            logS_list.append(logS)
            smiles_list.append(smiles)

    hist, edges, probs, bins_dict = get_bins_probs(logS_list, -11, 1, 1)
    for j in range(N):
        smiles_id = get_cmpd_idx(bins_dict)
        smiles = smiles_list[smiles_id]
        logS = logS_list[smiles_id]
        print(f"{smiles},{logS}")
