#! /usr/bin/env python

import sys
import matplotlib.pyplot as plt
import numpy as np


def read_measured_from_test(fname):
    data = {}
    with open(fname, "r") as fin:
        fin.readline()
        for line in fin:
            pairs = line.rstrip("\n").split(",")
            smiles = pairs[8]
            logS0 = float(pairs[2])
            data[smiles] = logS0

    return data


def read_predicted_values(fname):
    data = {}
    with open(fname, "r") as fin:
        for line in fin:
            pairs = line.rstrip("\n").split(",")
            smiles = pairs[0]
            logS0_pred = [float(x) for x in pairs[1:-1]]
            data[smiles] = logS0_pred

    return data


if __name__ == "__main__":
    smiles_pred = read_predicted_values(sys.argv[1])
    smiles_meas = read_measured_from_test(sys.argv[2])

    plot = False

    x = []
    y = []
    cols = []

    COUNT = 0
    RMSE = 0.0
    CORRECT = 0

    for key in smiles_meas.keys():
        pred_vals = smiles_pred[key]
        meas_val = smiles_meas[key]

        RMSE += (np.mean(pred_vals) - meas_val) ** 2.0
        if abs(np.mean(pred_vals) - meas_val) <= 0.5:
            CORRECT += 1
        COUNT += 1

        x.append(meas_val)
        y.append(np.mean(pred_vals))

    if plot:
        plt.scatter(x, y)
        plt.plot([-7, 0], [-7, 0], "--", color="grey", lw=2)
        plt.show()

    RMSE /= COUNT
    RMSE = np.sqrt(RMSE)

    CORRECT /= COUNT

    print(f"RMSE={RMSE}; \tCORRECT={CORRECT*100}%")
