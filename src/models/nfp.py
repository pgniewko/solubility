#!/usr/bin/env python
# Auhor: Pawel Gniewek, 2019-

import sys
import logging


from autograd import grad
import autograd.numpy as np
import autograd.numpy.random as npr

from neuralfingerprint import build_conv_deep_net
from neuralfingerprint import normalize_array, adam
from neuralfingerprint import build_batched_grad
from neuralfingerprint.util import rmse

from predictor import Predictor
from model_utils import get_training_data
from model_utils import parse_args

class NfpPredictor(Predictor):
    """
    TODO:
    """
    def __init__(self, radius=4, fplength=50):
        super().__init__()
        self._name = "NfpRegressor"

        self.model_params = dict(
                    fp_length=fplength,    # Usually neural fps need far fewer dimensions than morgan.
                    fp_depth=radius,       # The depth of the network equals the fingerprint radius.
                    conv_width=20,         # Only the neural fps need this parameter.
                    h1_size=100,            # Size of hidden layer of network on top of fps.
                    L2_reg=np.exp(-2))

        self.train_params = dict(
                    num_iters=100,
                    batch_size=100,
                    init_scale=np.exp(-4),
                    step_size=np.exp(-6))

        self.model = (None, None, None)

    def fit(self, smiles_list, logS_list, seed=0):
        train_smiles = [smiles for smiles in smiles_list]
        train_logS = [logS for logS in logS_list]

        conv_layer_sizes = [self.model_params['conv_width']] * self.model_params['fp_depth']
        conv_arch_params = {'num_hidden_features': conv_layer_sizes, 'fp_length': self.model_params['fp_length'], 'normalize': 1}

        # Neural net architecture
        net_arch_params = dict(layer_sizes = [self.model_params['fp_length'], self.model_params['h1_size']],
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
    args = parse_args()
    train_file = args.input
    results_file = args.output

    logging.basicConfig(stream=sys.stdout, level=logging.DEBUG)
    smiles_list, logS_list = get_training_data(train_file)
    nfp_regression = NfpPredictor()
    print(nfp_regression.train(smiles_list, logS_list, fname=results_file, y_randomization=args.y_rand))
    nfp_regression.plot()
