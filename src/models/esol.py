#!/usr/bin/env python
# Auhor: Pawel Gniewek, 2019-

import sys
import logging

from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen, Lipinski
from sklearn.linear_model import LinearRegression
from sklearn.linear_model import HuberRegressor
from sklearn.linear_model import RANSACRegressor
from sklearn.linear_model import TheilSenRegressor

from predictor import Predictor
from model_utils import get_training_data
from model_utils import parse_args

class ESOLCalculator(Predictor):
    """
    TODO:
    """

    def __init__(self, model='linear'):
        super().__init__()
        self._name = "ESOLCalculator"
        self.aromatic_query = Chem.MolFromSmarts("a")
        self._coef = {"MW": 0.0, "LogP": 0.0, "RB": 0.0, "AP": 0.0}
        self._intercept = 0.0

        # Reference params
        self._coef_esol = {"logP": -0.63, "MW": -0.0062, "RB": 0.066, "AP": -0.74}
        self._intercept_esol = 0.16
        self._coef_pat = {"logP": -0.74, "MW": -0.0066, "RB": 0.0034, "AP": -0.42}
        self._intercept_pat = 0.26

        self.model = model

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

    def fit(self, smiles_list, logS_list):
        X = []
        y = []
        for i, smiles in enumerate(smiles_list):
            mol = Chem.MolFromSmiles(smiles)
            (mw, logp, rotors, ap) = self._calc_esol_descriptors(mol)
            X.append([mw, logp, rotors, ap])
            y.append(logS_list[i])

        if self.model == 'linear':
            logging.debug(f'Model: {self.model}')
            model = LinearRegression()
        elif self.model == 'huber':
            logging.debug(f'Model: {self.model}')
            model = HuberRegressor(epsilon=1.5, alpha=2.0)
        elif self.model == 'ts':
            logging.debug(f'Model: {self.model}')
            model = TheilSenRegressor()
        else:
            logging.debug(f'Model: linear')
            model = LinearRegression()

        model.fit(X, y)
        self._intercept = model.intercept_
        self._coef["MW"] = model.coef_[0]
        self._coef["LogP"] = model.coef_[1]
        self._coef["RB"] = model.coef_[2]
        self._coef["AP"] = model.coef_[3]

    def predict(self, smiles_list):
        ypred = []
        for smiles in smiles_list:
            mol = Chem.MolFromSmiles(smiles)
            (mw, logp, rotors, ap) = self._calc_esol_descriptors(mol)
            y_val = self._intercept
            y_val += self._coef["MW"] * mw
            y_val += self._coef["LogP"] * logp
            y_val += self._coef["RB"] * rotors
            y_val += self._coef["AP"] * ap
            ypred.append(y_val)

        return ypred


if __name__ == "__main__":
    args = parse_args()
    train_file = args.input
    results_file = args.output

    logging.basicConfig(stream=sys.stdout, level=logging.DEBUG)
    smiles_list, logS_list = get_training_data(train_file)
    esol_calculator = ESOLCalculator(model=args.model)
    print(esol_calculator.train(smiles_list, logS_list, fname=results_file, y_randomization=args.y_rand))
    esol_calculator.plot(out_file=args.predictions_file)
