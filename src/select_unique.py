
import sys
import logging
from functools import partial
from collections import defaultdict

PROCESSED_PATH="/Users/pawel/Projects/solubility/data/processed"
TRAINING_DIR="/Users/pawel/Projects/solubility/data/training"



def one_outlier(num_list):
    tot = len(num_list)
    cnt = 0
    for el in num_list:
        if el > 0.0:
            cnt += 1

    if (tot - cnt) == 1 and tot > 1:
        return (True, 0)
    if cnt == 1 and tot > 1:
        return (True, 1)

    return (False,None)

def main():
    FILES_LIST = ['BOM.2017.JC.smi',
                  'D.2008.JCIC.smi',
                  'H.2000.JCIC.test1.smi',
                  'H.2000.JCIC.test2.smi',
                  'H.2000.JCIC.train.smi',
                  'HXZ.2004.JCIC.data_set.smi',
                  'HXZ.2004.JCIC.test.smi',
                  'LGG.2008.JCIM.100.smi',
                  'LGG.2008.JCIM.32.smi',
                  'OCHEM.Water.Solublity.05.27.2019.smi',
                  'POG.2007.JCIM.test.smi',
                  'POG.2007.JCIM.train.smi',
                  'WKH.2007.JCIM.smi']

    smiles_logS = defaultdict(list)
    for fname in FILES_LIST:
        logging.info(f"Processing {fname} file.")
        with open(f"{PROCESSED_PATH}/{fname}", 'r') as fin:
            for line in fin:
                pairs = line.rstrip('\n').split(',')
                smiles = pairs[0]
                logS = float(pairs[1])
                smiles_logS[smiles].append(logS)

    fout = f"{TRAINING_DIR}/solubility.uniq.smi"
    with open(fout, 'w', encoding='ascii') as fo:
        for key, values in smiles_logS.items():
            values.sort()
            vmin = values[0]
            vmax = values[-1]

            if vmin * vmax < 0:
                if len(values) == 2:
                    continue
            #    flag, minmax = one_outlier(values)
            #    if flag:
            #        if minmax == 0:
            #            values = values[1:]
            #        else:
            #            values = values[0:-1]
            
            if (vmax - vmin) > 1.0:
                continue

            fo.write(f"{key};")
            for val in values:
                fo.write(f"{val},")
            fo.write('\n')
                


if __name__ == "__main__":
    logging.basicConfig(stream=sys.stdout, level=logging.DEBUG)
    main()

 
