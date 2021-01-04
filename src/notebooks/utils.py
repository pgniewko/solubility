""" Utility functions, that for the purpose of this project,
are borrowed from Pat Walters: https://github.com/PatWalters/metk """

import math
import numpy as np
from scipy.stats import pearsonr
from scipy.stats import norm

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs
from rdkit.ML.Cluster import Butina


def pearson_confidence(r, num, interval=0.95):
    """
    Calculate upper and lower 95% CI for a Pearson r (not R**2)
    Inspired by https://stats.stackexchange.com/questions/18887
    :param r: Pearson's R
    :param num: number of data points
    :param interval: confidence interval (0-1.0)
    :return: lower bound, upper bound
    """
    stderr = 1.0 / math.sqrt(num - 3)
    z_score = norm.ppf(interval + (1.0 - interval) / 2.0)  # There was a bug in the original code
    delta = z_score * stderr
    lower = math.tanh(math.atanh(r) - delta)
    upper = math.tanh(math.atanh(r) + delta)
    return lower, upper


def max_possible_correlation(vals, error=0.6, method=pearsonr, cycles=1000):
    """
    Calculate the maximum possible correlation given a particular experimental error
    Based on Brown, Muchmore, Hajduk http://www.sciencedirect.com/science/article/pii/S1359644609000403
    :param vals: experimental values (should be on a log scale)
    :param error: experimental error
    :param method: method for calculating the correlation, must take 2 lists and return correlation and p_value
    :param cycles: number of random cycles
    :return: maximum possible correlation
    """
    cor_list = []
    for i in range(0, cycles):
        noisy_vals = []
        for val in vals:
            noisy_vals.append(val + np.random.normal(0, error))
        cor_list.append(method(vals, noisy_vals)[0])
    return np.mean(cor_list), np.array(cor_list)


def get_statisicts(statistic, models_data_files, statistics_map, RESULTS_PATH):
    """Read the statistics for each model
    :argument:
      statistic -- (string) models statistic
    :return:
      results -- (dict) a dictionary with statistic values (and standard deviations) for all available models.
    """
    results = {}
    for filename in models_data_files:
        with open('{}/{}'.format(RESULTS_PATH, filename), 'r') as fin:
            for line in fin:
                if line.startswith('#'):
                    continue
                pairs = line.rstrip('\n').split()

                model_idx = statistics_map[statistic]
                try:
                    results[pairs[0]] = [float(pairs[2*model_idx+1]), float(pairs[2*model_idx+2])]
                except ValueError as e:
                    print(e)
                    print("skipping the line")

    return results


def cluster_fps(fps, cutoff=0.2):
    """ Given the set of fingerpring return set of clusters.
    Similarity is calculated with Tanimoto metric.
    :arguments:
      fps (list) -- a list of fingerprints
      cutoff (flat) -- a cutoff for clustering. Compounds in the same claster have diff in similarity less than the cutoff

    :return:
      cs (list) -- list of clusters. Each cluster is given by a tuple with compounds ids as they occure in fps. The firs id is for the cluster exemplars.
    """
    # first generate the distance matrix:
    dists = []
    nfps = len(fps)
    for i in range(1, nfps):
        sims = DataStructs.BulkTanimotoSimilarity(fps[i], fps[:i])
        dists.extend([1 - x for x in sims])

    # now cluster the data:
    cs = Butina.ClusterData(dists, nfps, cutoff, isDistData=True)
    return cs


def get_data(fname, col=1, header=False):
    vals = []
    with open(fname, 'r') as fin:
        if header:
            fin.readline()
        for line in fin:
            pairs = line.rstrip('\n').split(',')
            vals.append(float(pairs[col]))

    return vals


def get_sim_data(ref_file, target_file):
    ref_smiles = []
    with open(ref_file, 'r') as fin:
        for line in fin:
            pairs = line.rstrip('\n').split(',')
            ref_smiles.append(pairs[0])

    target_smiles = []
    with open(target_file, 'r') as fin:
        for line in fin:
            pairs = line.rstrip('\n').split(',')
            target_smiles.append(pairs[0])

    ref_fps = [AllChem.GetMorganFingerprintAsBitVect(Chem.MolFromSmiles(x), 3, nBits=1024) for x in ref_smiles]
    target_fps = [AllChem.GetMorganFingerprintAsBitVect(Chem.MolFromSmiles(x), 3, nBits=1024) for x in target_smiles]

    mol_sims = []
    for i, fp_i in enumerate(target_fps):
        min_sim = 0.0
        for j, fp_j in enumerate(ref_fps):
            if i == j:
                continue
            sim_ij = DataStructs.TanimotoSimilarity(fp_i, fp_j)
            min_sim = max(min_sim, sim_ij)

        mol_sims.append(min_sim)

    return mol_sims
