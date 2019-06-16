#!/usr/bin/env python
# Auhor: Pawel Gniewek, 2019-

import sys
import logging

import numpy as np
from scipy import stats
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

from keras.models import Sequential
from keras.layers import Dense
from keras.layers import Dropout
from keras.regularizers import l2

from predictor import Predictor
from esol import ESOLCalculator
from rf import RFPredictor
from nfp import NfpPredictor


class EnsemblePredictor(Predictor):
    """
    """
    def __init__(self, fp='ecfp', radius=2, fp_length=1024, prop=False, n_ests=200):
        super().__init__()
        self._name = "EnsembleRegressor"
        self._fp = fp
        self._fp_r = radius
        self._fp_length = fp_length
        self._feats = list(rdMolDescriptors.Properties().GetPropertyNames())

        self.model = None
        self.esol_calculator = ESOLCalculator()
        self.rf_regression = RFPredictor(n_ests=n_ests, fp_type=self._fp)
        self.nfp_regression = NfpPredictor()

        self._means_logS = None
        self._std_logS = None
        self._means_props = None
        self._std_props = None
        self._size = -1
        self._epochs = 100


    def fit(self, smiles_list, logS_list):
        logging.info("Training ESOL model")
        self.esol_calculator.fit(smiles_list, logS_list)
        logging.info("Training RF model")
        self.rf_regression.fit(smiles_list, logS_list)
        logging.info("Training NFP model")
        self.nfp_regression.fit(smiles_list, logS_list)
       
        X_train = np.array(self._do_norm_X(smiles_list, find_norm=True))
        y_train = np.array(self._do_norm_y(logS_list))
        
        self.model = Sequential()
        self.model.add(Dense(50, 
                             input_dim=self._size, 
                             kernel_initializer='normal', 
                             activation='relu'))#, 
        self.model.add(Dense(25,
                             kernel_initializer='normal', 
                             activation='relu'))#,
        self.model.add(Dense(1, kernel_initializer='normal'))
        
        self.model.compile(loss='mean_squared_error', optimizer='adam')
        
        self.model.fit(X_train, y_train, epochs=self._epochs, batch_size=32,validation_split=0.0)

 
    def _do_norm_y(self, logS_list):
        y = []
        for i, smiles in enumerate(logS_list):
             y.append(logS_list[i])
        
        self._mean_logS = np.mean(y, axis=0)
        self._std_logS  = np.std(y, axis=0)
        
        logging.info("Means of Y")
        logging.info("Std of Y")

        y_norm = []
        for yi in y:
            yi_n = (yi - self._mean_logS) / self._std_logS
            y_norm.append ( yi_n)

        return y_norm
   
    def _undo_norm_y(self, y_pred):
        y = []
        for i, smiles in enumerate(y_pred):
             y.append(y_pred[i] * self._std_logS + self._mean_logS)
        return y

 
    def _do_norm_X(self, smiles_list, find_norm=True):
        X = []

        logS_list_esol = self.esol_calculator.predict(smiles_list)
        logS_list_rf   = self.rf_regression.predict(smiles_list)
        logS_list_nfp  = self.nfp_regression.predict(smiles_list)

        for i,smiles in enumerate(smiles_list):
            mol = Chem.MolFromSmiles(smiles)
            props = [list(rdMolDescriptors.Properties([name]).ComputeProperties(mol))[0] for name in self._feats] 
            vals = [logS_list_esol[i], logS_list_rf[i], logS_list_nfp[i]]
#            vals = [logS_list_esol[i], logS_list_nfp[i]]
            
            x_row = vals + props
            X.append(x_row)

        if find_norm:
            self._size = len(X[0])
            self._mean_props = np.mean(X, axis=0)
            self._std_props  = np.std(X, axis=0)
        
        logging.info("Size of props")
        print(self._size)
        logging.info("Means of X")
        print(self._mean_props)
        logging.info("Std of X")
        print(self._std_props)


        for i in range(len(X)):
            for j, X_ij in enumerate(X[i]):
                X[i][j] = (X_ij - self._mean_props[j]) / self._std_props[j]

        return X

    def _pickle(self, path, cv):
        pass

    def predict(self, smiles_list):     
        X = np.array( self._do_norm_X(smiles_list, find_norm=False))
        y_vals = self.model.predict(X)
        ypred = np.array(self._undo_norm_y(y_vals))
        return ypred

    
if __name__ == "__main__":
    from model_utils import get_training_data

    logging.basicConfig(stream=sys.stdout, level=logging.DEBUG)
    smiles_list, logS_list = get_training_data(sys.argv[1])
    enseble_regression = EnsemblePredictor(fp='maccs')
    print ( enseble_regression.train(smiles_list, logS_list) )
    enseble_regression.plot()

