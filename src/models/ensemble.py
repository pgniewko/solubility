#!/usr/bin/env python
# Auhor: Pawel Gniewek, 2019-


import sys
import logging

import numpy as np

from predictor import Predictor
from esol import ESOLCalculator
from rf import RFPredictor
from nfp import NfpPredictor


class EnsemblePredictor(Predictor):
    """
    """
    def __init__(self, fp='ecfp', radius=2, fp_length=1024, prop=False, n_ests=500):
        super().__init__()
        self._name = "EnsembleRegressor"
        self._fp = fp
        self._fp_r = radius
        self._fp_length = fp_length
        self._n_estimators = n_ests
        self.model = None
        self._prop = prop


    def fit(self, smiles_list, logS_list):
        pass


    def smiles_to_fps(self, smiles_list, fp_radius, fp_length):
        pass

# from rdkit.Chem import rdMolDescriptors
# self._feats = list(rdMolDescriptors.Properties().GetPropertyNames())
##:props = [list(rdMolDescriptors.Properties([name]).ComputeProperties(molecule))[0] for name in self._feats]

    def predict(self, smiles_list):
        pass

    
if __name__ == "__main__":
    from model_utils import get_training_data

    logging.basicConfig(stream=sys.stdout, level=logging.DEBUG)
    smiles_list, logS_list = get_training_data(sys.argv[1])
    enseble_regression = EnsemblePredictor()
    print ( enseble_regression.train(smiles_list, logS_list) )
    enseble_regression.plot()

