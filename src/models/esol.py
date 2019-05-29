#!/usr/bin/env python
# Auhor: Pawel Gniewek, 2019
# Code is based on the original Pat Walter's blog post:
# http://practicalcheminformatics.blogspot.com/2018/09/predicting-aqueous-solubility-its.html

import sys
from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen, Lipinski
from rdkit.Chem import PandasTools
import pandas as pd
from sklearn.linear_model import LinearRegression
from sklearn.model_selection import KFold
from scipy.stats import pearsonr
from sklearn.metrics import mean_squared_error
from sklearn.metrics import mean_absolute_error
from sklearn.metrics import r2_score
from collections import namedtuple
import numpy as np


class ESOLCalculator:
    """
    TODO: Write doc.
    """

    def __init__(self):
        self.aromatic_query = Chem.MolFromSmarts("a")
        self._coef = {"MW":0.0, "LogP":0.0, "RB":0.0, "AP":0.0}
        self._intercept = 0.0
        
        # As a reference
        self._coef_esol = {"logP": -0.63, "MW": -0.0062, "RB": 0.066, "AP": -0.74}
        self._intercept_esol = 0.16
        self._coef_pat  = {"logP": -0.74, "MW": -0.0066, "RB": 0.0034, "AP": -0.42}
        self._intercept_pat = 0.26


    def fit(self, smiles_list, logS_list):
        X = []
        y = []
        for i, smiles in enumerate(smiles_list):
            mol = Chem.MolFromSmiles(smiles)
            (mw, logp, rotors, ap) = self._calc_esol_descriptors(mol)
            X.append( [mw, logp, rotors, ap] )
            y.append(logS_list[i])

        model = LinearRegression()
        model.fit(X, y)
        coefficient_dict = {}
        self._intercept    = model.intercept_
        self._coef["MW"]   = model.coef_[0]
        self._coef["LogP"] = model.coef_[1]
        self._coef["RB"]   = model.coef_[2]
        self._coef["AP"]   = model.coef_[3]
    
    
    def predict(self, smiles_list):
        ypred = []
        for smiles in smiles_list:
            mol = Chem.MolFromSmiles(smiles)
            (mw, logp, rotors, ap) = self._calc_esol_descriptors(mol)
            y_val  = self._intercept
            y_val += self._coef["MW"] * mw
            y_val += self._coef["LogP"] * logp
            y_val += self._coef["RB"] * rotors
            y_val += self._coef["AP"] * ap 
            ypred.append(y_val)

        return ypred

    def score(self, smiles_list, logS_list, cv=5):
        scores = []
        
        kf = KFold(n_splits=cv)
        for train_index, test_index in kf.split(smiles_list):
            X_train = [smiles_list[idx] for idx in list(train_index)]
            y_train = [logS_list[idx] for idx in list(train_index)]
            X_test = [smiles_list[idx] for idx in list(test_index)]
            y_test = [logS_list[idx] for idx in list(test_index)]
            
            #X_train, X_test = smiles_list[train_index], smiles_list[test_index]
            #y_train, y_test = losS_list[train_index], logS_list[test_index]
            self.fit(X_train, y_train)
            scores.append( self.calc_score(X_test, y_test)  )


        return np.mean(scores, axis=0), np.std(scores, axis=0)


    def calc_score(self, smiles_list, y_true):
        y_pred = self.predict(smiles_list)
        mse = mean_squared_error(y_true, y_pred)
        mae = mean_absolute_error(y_true, y_pred)
        r2 = r2_score(y_true, y_pred)
        
        return (mse, mae, r2)


    def _calc_ap(self, mol):
        """
        Calculate aromatic proportion #aromatic atoms/#atoms total
        :param mol: input molecule
        :return: aromatic proportion
        """
        matches = mol.GetSubstructMatches(self.aromatic_query)
        return len(matches) / mol.GetNumAtoms()


    def _calc_esol_descriptors(self, mol):
        """
        Calcuate mw,logp,rotors and aromatic proportion (ap)
        :param mol: input molecule
        :return: tuple with calculated values
        """
        mw = Descriptors.MolWt(mol)
        logp = Crippen.MolLogP(mol)
        rotors = Lipinski.NumRotatableBonds(mol)
        ap = self._calc_ap(mol)
        return (mw, logp, rotors, ap)


def get_training_data(fname):
    cmpds_list = []
    logS_list = []
    with open(fname, 'r') as fin:
        for line in fin:
            pairs = line.rstrip('\n').split(',')
            smiles = pairs[0]
            logS = float(pairs[1])
            cmpds_list.append(smiles)
            logS_list.append(logS)

    return cmpds_list, logS_list


if __name__ == "__main__":

    smiles_list, logS_list = get_training_data(sys.argv[1])
    esol_calculator = ESOLCalculator()
    print ( esol_calculator.score(smiles_list, logS_list) )
    

#    esol_calculator.fit(smiles_list, logS_list)
#    y_pred = esol_calculator.predict(smiles_list)
#    esol_calculator.calc_score(smiles_list, logS_list)
#    import matplotlib.pyplot as plt
#    plt.plot(logS_list, y_pred, 'o', alpha=0.01)
#    plt.plot([-10,2],[-10,2], '--')
#    plt.show()

