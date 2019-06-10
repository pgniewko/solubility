import os
import logging
from abc import ABC, abstractmethod

import numpy as np
from sklearn.metrics import mean_squared_error
from sklearn.metrics import mean_absolute_error
from sklearn.metrics import r2_score
from sklearn.model_selection import KFold
import matplotlib.pyplot as plt


class Predictor(ABC):
    """
    """

    def __init__(self):
        self._logS_pred_data = []
        self._logS_exp_data = []
        pass

    @abstractmethod
    def fit(self):
        pass

    @abstractmethod
    def score(self):
        pass
    
    @abstractmethod
    def predict(self):
        pass
    
    @abstractmethod
    def _pickle(self, path, cv):
        pass

    def _make_dirs(self):
        return


    def train(self, smiles_list, logS_list, cv=5, save_flag=True, ensemble=False):
        run_path = os.path.dirname(os.path.abspath(__file__))
        saves_dir = f'{run_path}/saves'
        try:
            os.makedirs(saves_dir)
        except FileExistsError:
            logging.warning(f'Directory {saves_dir} exists')
 
        scores = []

        kf = KFold(n_splits=cv, shuffle=True, random_state=None)
        fold = 0
        for train_index, test_index in kf.split(smiles_list):
            logging.info('*{}* model is training {} fold'.format(self._name, fold))
            X_train = [smiles_list[idx] for idx in list(train_index)]
            y_train = [logS_list[idx] for idx in list(train_index)]
            X_test = [smiles_list[idx] for idx in list(test_index)]
            y_test = [logS_list[idx] for idx in list(test_index)]
            self.fit(X_train, y_train)
            self._pickle(saves_dir, fold)
            scores.append( self.score(X_test, y_test, save_flag)  )
            fold += 1

        return np.mean(scores, axis=0), np.std(scores, axis=0)


    def score(self, smiles_list, y_true, save_flag=True):
        y_pred = self.predict(smiles_list)
        mse = mean_squared_error(y_true, y_pred)
        mae = mean_absolute_error(y_true, y_pred)
        r2 = r2_score(y_true, y_pred)
        
        if save_flag:
            self._logS_pred_data += list(y_pred)
            self._logS_exp_data  += list(y_true)

        return (mse, mae, r2)

    def plot(self):
        plt.figure(figsize=(7,7))
        plt.plot(self._logS_exp_data, self._logS_pred_data, 'o',alpha=0.05)
        plt.xlabel('logS0 [mol/L] (measured)', fontsize=14)
        plt.ylabel('logS0 [mol/L] (predicted)', fontsize=14)
        plt.title('Model: {}'.format(self._name))
        plt.show()
