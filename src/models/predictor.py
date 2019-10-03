import logging
from abc import ABC, abstractmethod

import numpy as np
from sklearn.metrics import mean_squared_error
from sklearn.metrics import mean_absolute_error
from sklearn.metrics import r2_score
from sklearn.model_selection import KFold
import matplotlib.pyplot as plt

from model_utils import bootstrap_error_estimate

class Predictor(ABC):
    """
    TODO:
    """

    def __init__(self):
        self._logS_pred_data = []
        self._logS_exp_data = []
        self._name = None
        pass

    @abstractmethod
    def fit(self):
        pass

    @abstractmethod
    def predict(self):
        pass

    def get_name(self):
        return self._name

    def train(self, train_smiles, logS_list, cv=5, fname=None):
        scores = []

        kf = KFold(n_splits=cv, shuffle=True, random_state=None)
        fold = 0
        for train_index, validate_index in kf.split(train_smiles):
            logging.info('*{}* model is training {} fold'.format(self._name, fold))
            X_train = [train_smiles[idx] for idx in list(train_index)]
            y_train = [logS_list[idx] for idx in list(train_index)]
            X_validate = [train_smiles[idx] for idx in list(validate_index)]
            y_validate = [logS_list[idx] for idx in list(validate_index)]
            self.fit(X_train, y_train)
            scores.append(self.score(X_validate, y_validate))
            fold += 1

        if fname != None:
            with open(fname, 'w') as fout:
                fout.write(f'{self._name}\t')
                means = np.mean(scores, axis=0)
                stds = np.std(scores, axis=0)
                for i, mean_ in enumerate(means):
                    fout.write('{}\t{}\t'.format(mean_, stds[i]))

                mse_glob = mean_squared_error(self._logS_pred_data, self._logS_pred_data)
                mae_glob = mean_absolute_error(self._logS_pred_data, self._logS_pred_data)
                r2_glob = r2_score(self._logS_pred_data, self._logS_pred_data)
                mse_lower, mse_upper = bootstrap_error_estimate(self._logS_pred_data,
                                                                self._logS_pred_data,
                                                                method=mean_squared_error)
                mae_lower, mae_upper = bootstrap_error_estimate(self._logS_pred_data,
                                                                self._logS_pred_data,
                                                                method=mean_absolute_error)
                
                r2_lower, r2_upper = bootstrap_error_estimate(self._logS_pred_data,
                                                                self._logS_pred_data,
                                                                method=r2_score)
                four.write(f'{mse_glob}\t{mae_lower}\t{mse_upper}\t')
                four.write(f'{mae_glob}\t{mae_lower}\t{mae_upper}\t')
                four.write(f'{r2_glob}\t{r2_lower}\t{r2_upper}\t')
                fout.write('\n')

        return np.mean(scores, axis=0), np.std(scores, axis=0)


    def test(self, train_smiles, logS_list, test_smiles, cv=5):
        predictions = []

        kf = KFold(n_splits=cv, shuffle=True, random_state=None)
        fold = 0
        for train_index, validate_index in kf.split(train_smiles):
            logging.info('*{}* model is training {} fold'.format(self._name, fold))
            X_train = [train_smiles[idx] for idx in list(train_index)]
            y_train = [logS_list[idx] for idx in list(train_index)]
#            X_validate = [train_smiles[idx] for idx in list(validate_index)]
#            y_validate = [logS_list[idx] for idx in list(validate_index)]
            self.fit(X_train, y_train)

            # Don't do validation - go to testing with the model at hand
            y_pred = self.predict(test_smiles)
            predictions.append(y_pred)
            fold += 1

        return predictions

    def score(self, smiles_list, y_true, save_flag=True):
        y_pred = self.predict(smiles_list)
        mse = mean_squared_error(y_true, y_pred)
        mae = mean_absolute_error(y_true, y_pred)
        r2 = r2_score(y_true, y_pred)

        if save_flag:
            self._logS_pred_data += list(y_pred)
            self._logS_exp_data += list(y_true)

        return (mse, mae, r2)

    def plot(self):
        plt.figure(figsize=(7, 7))
        plt.plot(self._logS_exp_data, self._logS_pred_data, 'o', alpha=0.05)
        plt.xlabel('logS0 [mol/L] (measured)', fontsize=14)
        plt.ylabel('logS0 [mol/L] (predicted)', fontsize=14)
        plt.title('Model: {}'.format(self._name))
        plt.show()
