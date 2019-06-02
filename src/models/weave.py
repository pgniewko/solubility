#!/usr/bin/env python
# Auhor: Pawel Gniewek, 2019-


import sys
import logging

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
import deepchem
from deepchem.models import WeaveModel
from rdkit import Chem, DataStructs
from sklearn import preprocessing

from predictor import Predictor


class WeavePredictor(Predictor):
    """
    """

    def __init__(self):
        self._name = "WeaveRegressor"
        self._featurizer = deepchem.feat.WeaveFeaturizer()
        self.model = None


    def fit(self, smiles_list, logS_list):
        featurizer = self._featurizer #deepchem.feat.WeaveFeaturizer()

        d = {'SMILES': smiles_list, 'logS': logS_list}
        df = pd.DataFrame(data=d)
        dataset = featurize_smiles_df(df, featurizer, "SMILES", log_every_N=1000, verbose=True)                
 
        move_mean = True
        transformers = [deepchem.trans.NormalizationTransformer(transform_y=True, dataset=dataset, move_mean=True)]
        
        for transformer in transformers:
             dataset = transformer.transform(dataset)
             
        return dataset, featurizer, transformers


def run_model(model_func, task_list, featurizer, normalize, training_file_name, validation_file_name, nb_epoch):
    dataset, featurizer, transformers = featurize_data(task_list, featurizer, normalize, training_file_name)
    model = model_func()
    if nb_epoch > 0:
        model.fit(dataset, nb_epoch)
    else:
        model.fit(dataset)
    pred_df = generate_prediction(validation_file_name, model, featurizer, transformers)
    return pred_df


        X = self.smiles_to_fps(smiles_list, 3, self._fp_length)
        y = [logS for logS in logS_list]


        batch_size = 64
        self.model = WeaveModel(1, batch_size=batch_size, learning_rate=1e-3, use_queue=False, mode='regression')
       
        nb_epoch = 30
        self.model.fit(dataset, nb_epoch)
        
    def featurize_data(tasks, featurizer, normalize, dataset_file):
        loader = deepchem.data.CSVLoader(tasks=tasks, smiles_field="SMILES", featurizer=featurizer)
        dataset = loader.featurize(dataset_file, shard_size=8192)
        move_mean = True
        if normalize:
            transformers = [deepchem.trans.NormalizationTransformer(
                transform_y=True, dataset=dataset, move_mean=move_mean)]
        else:
            transformers = []
        for transformer in transformers:
             dataset = transformer.transform(dataset)
        
        return dataset, featurizer, transformers
    
    
    
    def predict(self, smiles_list):
        X = self.smiles_to_fps(smiles_list, 2, self._fp_length)
        #X = self._scaler.transform(X)
        y_pred = self.model.predict(X)
        return y_pred

    
if __name__ == "__main__":
    from model_utils import get_training_data

    logging.basicConfig(stream=sys.stdout, level=logging.DEBUG)
    smiles_list, logS_list = get_training_data(sys.argv[1])
    rf_regression = MlpPredictor()
    print ( rf_regression.train(smiles_list, logS_list) )
    
