#! /usr/bin/env python

import sys
import logging

from models.esol import ESOLCalculator
from models.rf import RFPredictor
from models.nfp import NfpPredictor
from models.ensemble import EnsemblePredictor
from models.model_utils import get_training_data
from models.model_utils import get_test_data

def save_predictions(test_smiles, results_data, fname):
    num_folds = len(results_data)
    with open(fname, 'w') as fout:
        for i, smiles in enumerate(test_smiles):
            fout.write(f'{smiles} ')
            for j in range(num_folds):
                fout.write("{} ".format(results_data[j][i]))
            fout.write('\n')
                

if __name__ == "__main__":
    
    logging.basicConfig(stream=sys.stdout, level=logging.DEBUG)
    train_smiles, logS_list = get_training_data(sys.argv[1])
    test_smiles = get_test_data(sys.argv[2])
    
    production_model = EnsemblePredictor(fp='maccs')
    results_data = production_model.test(train_smiles, logS_list, test_smiles)
    save_predictions(test_smiles, results_data, sys.argv[3], m_name)
   

