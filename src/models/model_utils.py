import argparse

def get_training_data(fname):
    """
    """

    cmpds_list = []
    logS_list = []
    with open(fname, 'r') as fin:
        for line in fin:
            pairs = line.rstrip('\n').split(',')
            smiles = pairs[0]
            logS = float(pairs[1])
            cmpds_list.append(smiles)
            logS_list.append(logS)

    return cmpds_list, logS_list


def get_test_data(fname):
    """
    """

    cmpds_list = []
    with open(fname, 'r') as fin:
        for line in fin:
            smiles = line.rstrip('\n')
            cmpds_list.append(smiles)

    return cmpds_list

def parse_args():
    parser = argparse.ArgumentParser(description='__doc__')
    parser.add_argument('--input', help='Training file')
    parser.add_argument('--output', help='Output file')
    parser.add_argument('--predictions_file', help='Save measurements and predictions')
    parser.add_argument('--model', help='Regression model: linear, huber [default: linear]', default=None)
    parser.add_argument('--y_rand', default=False, action='store_true', help='Run in y-randomization mode')
    args = parser.parse_args()
    return args

