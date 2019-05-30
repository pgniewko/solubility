import logging
from abc import ABC, abstractmethod

import numpy as np
from sklearn.metrics import mean_squared_error
from sklearn.metrics import mean_absolute_error
from sklearn.metrics import r2_score
from sklearn.model_selection import KFold


class Predictor(ABC):
    """
    """

    def __init__(self):
        self._name = None
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
    
    def _make_dirs(self):
        return

    def _pickle(self):
        return

    def train(self, smiles_list, logS_list, cv=5):
        scores = []

        kf = KFold(n_splits=cv, shuffle=True, random_state=None)
        fold = 0
        for train_index, test_index in kf.split(smiles_list):
            print
            logging.info('*{}* model is training {} fold'.format(self._name, fold))
            X_train = [smiles_list[idx] for idx in list(train_index)]
            y_train = [logS_list[idx] for idx in list(train_index)]
            X_test = [smiles_list[idx] for idx in list(test_index)]
            y_test = [logS_list[idx] for idx in list(test_index)]
            self.fit(X_train, y_train)
            scores.append( self.score(X_test, y_test)  )
            fold += 1

        return np.mean(scores, axis=0), np.std(scores, axis=0)


    def score(self, smiles_list, y_true):
        y_pred = self.predict(smiles_list)
        mse = mean_squared_error(y_true, y_pred)
        mae = mean_absolute_error(y_true, y_pred)
        r2 = r2_score(y_true, y_pred)

        return (mse, mae, r2)

        
