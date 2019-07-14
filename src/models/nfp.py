#!/usr/bin/env python
# Auhor: Pawel Gniewek, 2019-

import sys
import logging

from predictor import Predictor

from autograd import grad
import autograd.numpy as np
import autograd.numpy.random as npr

from neuralfingerprint import build_conv_deep_net
from neuralfingerprint import normalize_array, adam
from neuralfingerprint import build_batched_grad
from neuralfingerprint.util import rmse


class NfpPredictor(Predictor):
    """
    """
    def __init__(self, radius=4, fplength=64):
        super().__init__()
        self._name = "NfpRegressor"

        self.model_params = dict(
                    fp_length=fplength,    # Usually neural fps need far fewer dimensions than morgan.
                    fp_depth=radius,       # The depth of the network equals the fingerprint radius.
                    conv_width=20,         # Only the neural fps need this parameter.
                    h1_size=64,            # Size of hidden layer of network on top of fps.
#                    h2_size=128,           # Size of hidden layer of network on top of fps.
                    L2_reg=np.exp(-2))

        self.train_params = dict(
                    num_iters=50,
                    batch_size=64,
                    init_scale=np.exp(-4),
                    step_size=np.exp(-6))

        self.model = (None, None, None)

    def fit(self, smiles_list, logS_list, seed=0):
        train_smiles = [smiles for smiles in smiles_list]
        train_logS = [logS for logS in logS_list]

        conv_layer_sizes = [self.model_params['conv_width']] * self.model_params['fp_depth']
        conv_arch_params = {'num_hidden_features': conv_layer_sizes, 'fp_length': self.model_params['fp_length'], 'normalize': 1}

        # Neural net architecture
        net_arch_params = dict(
                layer_sizes = [
                         self.model_params['fp_length'],
                         self.model_params['h1_size']],
#                         self.model_params['h2_size']],  # Two hidden layers
                         normalize=True,
                         L2_reg = self.model_params['L2_reg'],
                         nll_func = rmse)

        loss_fun, pred_fun, conv_parser = build_conv_deep_net(conv_arch_params, net_arch_params, self.model_params['L2_reg'])

        num_weights = len(conv_parser)
        init_weights = npr.RandomState(seed).randn(num_weights) * self.train_params['init_scale']

        train_logS_norm, undo_norm = normalize_array(train_logS)

        # Build gradient using autograd.
        grad_fun = grad(loss_fun)
        grad_fun_with_data = build_batched_grad(
                grad_fun,
                self.train_params['batch_size'],
                train_smiles,
                train_logS_norm)

        # Optimize weights.
        trained_weights = adam(
                grad_fun_with_data,
                init_weights,
                num_iters=self.train_params['num_iters'],
                step_size=self.train_params['step_size'])

        self.model = (undo_norm, trained_weights, pred_fun)

    def predict(self, smiles_list):
        (undo_norm, trained_weights, pred_fun) = self.model
        y_pred = undo_norm(pred_fun(trained_weights, smiles_list))
        return y_pred


if __name__ == "__main__":
    from model_utils import get_training_data

    logging.basicConfig(stream=sys.stdout, level=logging.DEBUG)
    smiles_list, logS_list = get_training_data(sys.argv[1])
    nfp_regression = NfpPredictor()
    print(nfp_regression.train(smiles_list, logS_list))
    nfp_regression.plot()
