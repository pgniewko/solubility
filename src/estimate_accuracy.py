#! /usr/bin/env python

import sys
import matplotlib.pyplot as plt
import numpy as np


def read_smiles_and_gse(fname):
    data = {}
    with open(fname, 'r') as fin:
        for line in fin:
            pairs = line.rstrip('\n').split(',')
            smiles = pairs[0]
            gse = float(pairs[1])
            data[smiles] = gse

    return data


def read_measured_from_test(fname):
    data = {}
    with open(fname, 'r') as fin:
        for line in fin:
            pairs = line.rstrip('\n').split(',')
            smiles = pairs[0]
            logS0 = float(pairs[1])
            data[smiles] = logS0

    return data


def read_predicted_values(fname):
    data = {}
    with open(fname, 'r') as fin:
        for line in fin:
            pairs = line.rstrip('\n').split(',')
            smiles = pairs[0]
            print(pairs[1:-1])
            logS0_pred = [float(x) for x in pairs[1:-1]]
            data[smiles] = logS0_pred

    return data


if __name__ == "__main__":
    smiles_gse = read_smiles_and_gse(sys.argv[1])
    smiles_pred = read_predicted_values(sys.argv[2])
    smiles_meas = read_measured_from_test(sys.argv[3])

    x = []
    y = []
    cols = []

    cnt = 0
    rmsd_gse = 0.0
    rmsd_pg = 0.0

    for key in smiles_meas.keys():
        gse_val = smiles_gse[key]
        pred_vals = smiles_pred[key]
        meas_val = smiles_meas[key]

        rmsd_gse += (gse_val - meas_val)**2.0
        rmsd_pg += (np.mean(pred_vals) - meas_val)**2.0

        x.append(meas_val)
        y.append(gse_val)
        cols.append('red')

        for val_ in pred_vals:
            x.append(meas_val)
            y.append(val_)
            cols.append('blue')

        x.append(meas_val)
        y.append(np.mean(pred_vals))
        cols.append('green')
        cnt += 1

    plt.scatter(x, y, c=cols)
    plt.plot([-7, 0], [-7, 0], '--', color='grey', lw=2)
    plt.show()

    rmsd_gse /= cnt
    rmsd_gse = np.sqrt(rmsd_gse)
    rmsd_pg /= cnt
    rmsd_pg = np.sqrt(rmsd_pg)

    print("GSE (RMSD) = {}\nPG(RMSD) = {}\ncnt={}".format(rmsd_gse, rmsd_pg, cnt))
