def get_training_data(fname):
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
    cmpds_list = []
    with open(fname, 'r') as fin:
        for line in fin:
            smiles = line.rstrip('\n')
            cmpds_list.append(smiles)

    return cmpds_list
