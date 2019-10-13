#! /usr/bin/env python

import sys
import logging
import argparse

from models.rf import RFPredictor
from models.nfp import NfpPredictor
from models.esol import ESOLCalculator
from models.ensemble import EnsemblePredictor

from models.model_utils import get_training_data
from models.model_utils import get_test_data


def save_predictions(test_smiles, results_data, fname):
    """
    TODO: other statitics are going to be returned with results_data.
    results_data has to be a dict.
    """

    num_folds = len(results_data)
    with open(fname, 'w') as fout:
        for i, smiles in enumerate(test_smiles):
            fout.write(f'{smiles},')
            for j in range(num_folds):
                fout.write("{},".format(float(results_data[j][i])))
            fout.write('\n')


def main(args):
    logging.basicConfig(stream=sys.stdout, level=logging.DEBUG)
    train_smiles, train_logS = get_training_data(args.train_file)
    test_smiles = get_test_data(args.test_file)

    if args.model == 'esol':
        production_model = ESOLCalculator()
    elif args.model == 'rf':
        production_model = RFPredictor()
    elif args.model == 'nfp':
        production_model = NfpPredictor()
    elif args.model == 'ensemble':
        production_model = EnsemblePredictor()
    else:
        print("Model not recognized: {}".format(args.model))
        return 1

    results_data = production_model.test(train_smiles, train_logS, test_smiles)
    save_predictions(test_smiles, results_data, args.out_file)

    return 0


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--train_file', help='Path to the file with the train data')
    parser.add_argument('--test_file', help='Path to the file with the test data')
    parser.add_argument('--out_file', help='Path to the output file')
    parser.add_argument('--model', default='rf', help='Choose the model: [esol, rf, nfp, ensemble]')
    args = parser.parse_args()
    return args


if __name__ == "__main__":
    args = parse_args()
    sys.exit(main(args))
