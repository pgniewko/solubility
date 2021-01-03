import os
import sys
import logging
from collections import defaultdict

import numpy as np
from utils import canonicalize_smiles, mol_wt
from rdkit import Chem


# Get rid of the hard coded paths
DATA_PATH = "/Users/pawel/projects/solubility/data/raw"
PROCESSED_PATH = "/Users/pawel/projects/solubility/data/processed"
TRAINING_DIR = "/Users/pawel/projects/solubility/data/training"
TEST_PATH = "/Users/pawel/projects/solubility/data/test"


def process_AB_2001_EJPS():
    fname = f"{DATA_PATH}/AB.2001.EJPS.txt"
    fout = f"{PROCESSED_PATH}/AB.2001.EJPS.smi"
    logging.info(f"Processing {fname}")

    cmpd_list = []
    with open(fname, 'r') as fin:
        for line in fin:
            if not line.startswith("#"):
                pairs = line.rstrip('\n').split(",")
                canon_smiles = canonicalize_smiles(pairs[1])
                logS = float(pairs[2])
                cmpd_list.append((canon_smiles, logS))

    with open(fout, 'w') as fo:
        for el in cmpd_list:
            smiles = el[0]
            logS = el[1]
            fo.write(f"{smiles},{logS}\n")

    logging.info(f"Saved {fout}")


def process_ABB_2000_PR():
    fname = f"{DATA_PATH}/ABB.2000.PR.txt"
    fout = f"{PROCESSED_PATH}/ABB.2000.PR.smi"
    logging.info(f"Processing {fname}")

    cmpd_list = []
    with open(fname, 'r') as fin:
        for line in fin:
            if not line.startswith("#"):
                pairs = line.rstrip('\n').split(",")
                canon_smiles = canonicalize_smiles(pairs[1])
                logS = float(pairs[2])
                cmpd_list.append((canon_smiles, logS))

    with open(fout, 'w') as fo:
        for el in cmpd_list:
            smiles = el[0]
            logS = el[1]
            fo.write(f"{smiles},{logS}\n")

    logging.info(f"Saved {fout}")


def process_BOM_2017_JC():
    fname = f"{DATA_PATH}/BOM.2017.JC.txt"
    fout = f"{PROCESSED_PATH}/BOM.2017.JC.smi"
    logging.info(f"Processing {fname}")

    cmpd_list = []
    with open(fname, 'r') as fin:
        for line in fin:
            pairs = line.rstrip('\n').split(",")
            canon_smiles = canonicalize_smiles(pairs[0])
            logS = float(pairs[1])
            cmpd_list.append((canon_smiles, logS))

    with open(fout, 'w', encoding='ascii') as fo:
        for el in cmpd_list:
            smiles = el[0]
            logS = el[1]
            fo.write(f"{smiles},{logS}\n")

    logging.info(f"Saved {fout}")


def process_D_2008_JCIC():
    fname = f"{DATA_PATH}/D.2008.JCIC.solubility.v1.txt"
    fout = f"{PROCESSED_PATH}/D.2008.JCIC.smi"
    logging.info(f"Processing {fname}")

    cmpd_list = []
    with open(fname, 'r') as fin:
        for line in fin:
            if not line.startswith("#"):
                pairs = line.rstrip('\n').split(",")
                canon_smiles = canonicalize_smiles(pairs[-1])
                logS = float(pairs[-3])
                cmpd_list.append((canon_smiles, logS))

    with open(fout, 'w', encoding='ascii') as fo:
        for el in cmpd_list:
            smiles = el[0]
            logS = el[1]
            fo.write(f"{smiles},{logS}\n")

    logging.info(f"Saved {fout}")


def process_H_2000_JCIC_test1():
    fname = f"{DATA_PATH}/H.2000.JCIC.test1.txt"
    fout = f"{PROCESSED_PATH}/H.2000.JCIC.test1.smi"
    logging.info(f"Processing {fname}")

    cmpd_list = []
    with open(fname, 'r', encoding="ISO-8859-1") as fin:
        for line in fin:
            if not line.startswith("#"):
                try:
                    pairs = line.rstrip('\n').split()
                    canon_smiles = canonicalize_smiles(pairs[6])
                    logS = float(pairs[3])
                    cmpd_list.append((canon_smiles, logS))
                except:
                    smiles = pairs[6]
                    logging.info(f"Failed to process {smiles}")

    with open(fout, 'w', encoding='ascii') as fo:
        for el in cmpd_list:
            smiles = el[0]
            logS = el[1]
            fo.write(f"{smiles},{logS}\n")

    logging.info(f"Saved {fout}")


def process_H_2000_JCIC_test2():
    fname = f"{DATA_PATH}/H.2000.JCIC.test2.txt"
    fout = f"{PROCESSED_PATH}/H.2000.JCIC.test2.smi"
    logging.info(f"Processing {fname}")

    cmpd_list = []
    with open(fname, 'r', encoding="ISO-8859-1") as fin:
        for line in fin:
            if not line.startswith("#"):
                try:
                    pairs = line.rstrip('\n').split()
                    canon_smiles = canonicalize_smiles(pairs[6], SLN=True)
                    logS = float(pairs[3])
                    cmpd_list.append((canon_smiles, logS))
                except:
                    smiles = pairs[6]
                    logging.info(f"Failed to process {smiles}")

    with open(fout, 'w', encoding='ascii') as fo:
        for el in cmpd_list:
            smiles = el[0]
            logS = el[1]
            fo.write(f"{smiles},{logS}\n")

    logging.info(f"Saved {fout}")


def process_H_2000_JCIC_train():
    fname = f"{DATA_PATH}/H.2000.JCIC.train.txt"
    fout = f"{PROCESSED_PATH}/H.2000.JCIC.train.smi"
    logging.info(f"Processing {fname}")

    cmpd_list = []
    with open(fname, 'r', encoding="ISO-8859-1") as fin:
        for line in fin:
            if not line.startswith("#"):
                try:
                    pairs = line.rstrip('\n').split()
                    canon_smiles = canonicalize_smiles(pairs[6])
                    logS = float(pairs[3])
                    cmpd_list.append((canon_smiles, logS))
                except:
                    smiles = pairs[6]
                    logging.info(f"Failed to process {smiles}")

    with open(fout, 'w', encoding="ascii") as fo:
        for el in cmpd_list:
            smiles = el[0]
            logS = el[1]
            fo.write(f"{smiles},{logS}\n")

    logging.info(f"Saved {fout}")


def process_HXZ_2004_JCIC_data():
    fname = f"{DATA_PATH}/HXZ.2004.JCIC.data_set.txt"
    fout = f"{PROCESSED_PATH}/HXZ.2004.JCIC.data_set.smi"
    logging.info(f"Processing {fname}")

    cmpd_list = []
    with open(fname, 'r') as fin:
        for line in fin:
            if not line.startswith("#"):
                try:
                    pairs = line.rstrip('\n').split()
                    canon_smiles = canonicalize_smiles(pairs[0])
                    logS = float(pairs[2])
                    cmpd_list.append((canon_smiles, logS))
                except:
                    smiles = pairs[0]
                    logging.info(f"Failed to process {smiles}")

    with open(fout, 'w', encoding="ascii") as fo:
        for el in cmpd_list:
            smiles = el[0]
            logS = el[1]
            fo.write(f"{smiles},{logS}\n")

    logging.info(f"Saved {fout}")


def process_LGG_2008_JCIM_100():
    fname = f"{DATA_PATH}/LGG.2008.JCIM.100.txt"
    fout = f"{PROCESSED_PATH}/LGG.2008.JCIM.100.smi"
    logging.info(f"Processing {fname}")

    cmpd_list = []
    with open(fname, 'r') as fin:
        cnt = 0
        for line in fin:
            if not line.startswith("#"):
                pairs = line.rstrip('\n').split(',')
                canon_smiles = canonicalize_smiles(pairs[0])
                logS = float(pairs[2])
                cmpd_list.append((canon_smiles, logS))
                cnt += 1

    with open(fout, 'w', encoding="ascii") as fo:
        for el in cmpd_list:
            smiles = el[0]
            logS = el[1]
            fo.write(f"{smiles},{logS}\n")

    logging.info(f"Saved {fout}")


def process_LGG_2008_JCIM_32():
    fname = f"{DATA_PATH}/LGG.2008.JCIM.32.txt"
    fout = f"{PROCESSED_PATH}/LGG.2008.JCIM.32.smi"
    logging.info(f"Processing {fname}")

    cmpd_list = []
    with open(fname, 'r') as fin:
        cnt = 0
        for line in fin:
            if not line.startswith("#"):
                pairs = line.rstrip('\n').split(',')
                canon_smiles = canonicalize_smiles(pairs[0])
                mg = float(pairs[1]) / 1000
                try:
                    mw = mol_wt(canon_smiles)
                except TypeError as e:
                    logging.error(f"TypeError for {canon_smiles}: {e}")
                    continue
                logS = np.log10(mg/mw)
                cmpd_list.append((canon_smiles, logS))
                cnt += 1

    with open(fout, 'w', encoding="ascii") as fo:
        for el in cmpd_list:
            smiles = el[0]
            logS = el[1]
            fo.write(f"{smiles},{logS}\n")

    logging.info(f"Saved {fout}")


def process_POG_2007_JCIM_test():
    fname = f"{DATA_PATH}/POG.2007.JCIM.test.txt"
    fout = f"{PROCESSED_PATH}/POG.2007.JCIM.test.smi"
    logging.info(f"Processing {fname}")

    cmpd_list = []
    with open(fname, 'r') as fin:
        for line in fin:
            if not line.startswith("#"):
                try:
                    pairs = line.rstrip('\n').split(";")
                    canon_smiles = canonicalize_smiles(pairs[1].strip("\""))
                    logS = float(pairs[4])
                    cmpd_list.append((canon_smiles, logS))
                except:
                    smiles = pairs[1]
                    logging.info(f"Failed to process {smiles}")

    with open(fout, 'w', encoding="ascii") as fo:
        for el in cmpd_list:
            smiles = el[0]
            logS = el[1]
            fo.write(f"{smiles},{logS}\n")

    logging.info(f"Saved {fout}")


def process_POG_2007_JCIM_train():
    fname = f"{DATA_PATH}/POG.2007.JCIM.train.txt"
    fout = f"{PROCESSED_PATH}/POG.2007.JCIM.train.smi"
    logging.info(f"Processing {fname}")

    cmpd_list = []
    with open(fname, 'r') as fin:
        for line in fin:
            if not line.startswith("#"):
                try:
                    pairs = line.rstrip('\n').split(";")
                    canon_smiles = canonicalize_smiles(pairs[1].strip("\""))
                    logS = float(pairs[4])
                    cmpd_list.append((canon_smiles, logS))
                except:
                    smiles = pairs[1]
                    logging.info(f"Failed to process {smiles}")

    with open(fout, 'w', encoding="ascii") as fo:
        for el in cmpd_list:
            smiles = el[0]
            logS = el[1]
            fo.write(f"{smiles},{logS}\n")

    logging.info(f"Saved {fout}")


def process_HXZ_2004_JCIC_test():
    fname = f"{DATA_PATH}/HXZ.2004.JCIC.test_set1.sdf"
    fout = f"{PROCESSED_PATH}/HXZ.2004.JCIC.test.smi"
    logging.info(f"Processing {fname}")

    suppl = Chem.SDMolSupplier(fname)
    with open(fout, 'w', encoding='ascii') as fo:
        for mol in suppl:
            smiles = canonicalize_smiles(Chem.MolToSmiles(mol))
            logS = str(mol.GetProp('logS'))
            fo.write(f'{smiles},{logS}\n')

    logging.info(f"Saved to {fout}")


def process_WKH_2007_JCIM():
    fname = f"{DATA_PATH}/WKH.2007.JCIM.solubility.sdf"
    fout = f"{PROCESSED_PATH}/WKH.2007.JCIM.smi"
    logging.info(f"Processing {fname}")

    suppl = Chem.SDMolSupplier(fname)
    with open(fout, 'w', encoding='ascii') as fo:
        for mol in suppl:
            smiles = canonicalize_smiles(Chem.MolToSmiles(mol))
            logS = str(mol.GetProp('EXPT'))
            fo.write(f'{smiles},{logS}\n')

    logging.info(f"Saved to {fout}")


def process_WHX_2009_JCIM():
    files_in = ["WHX.2009.JCIM.Set-001.csv",
                "WHX.2009.JCIM.Set-002.csv",
                "WHX.2009.JCIM.Set-003.csv",
                "WHX.2009.JCIM.Set-004.csv",
                "WHX.2009.JCIM.Set-005.csv"]

    files_out = ["WHX.2009.JCIM.Set-001.smi",
                 "WHX.2009.JCIM.Set-002.smi",
                 "WHX.2009.JCIM.Set-003.smi",
                 "WHX.2009.JCIM.Set-004.smi",
                 "WHX.2009.JCIM.Set-005.smi"]

    for i, file_in in enumerate(files_in):
        fname = f"{DATA_PATH}/{file_in}"
        file_out = files_out[i]
        fout = f"{PROCESSED_PATH}/{file_out}"
        logging.info(f"Processing {fname}")

        cmpd_list = []
        with open(fname, 'r') as fin:
            for line in fin:
                if not line.startswith("#"):
                    try:
                        pairs = line.rstrip('\n').split(",")
                        canon_smiles = canonicalize_smiles(pairs[1], SLN=True)
                        logS = float(pairs[0])
                        cmpd_list.append((canon_smiles, logS))
                    except:
                        smiles = pairs[1]
                        logging.info(f"Failed to process {smiles}")

        with open(fout, 'w', encoding="ascii") as fo:
            for el in cmpd_list:
                smiles = el[0]
                logS = el[1]
                fo.write(f"{smiles},{logS}\n")

        logging.info(f"Saved {fout}")


def process_OCHEM():
    fname = f"{DATA_PATH}/OCHEM.Water.Solublity.05.27.2019.txt"
    fout = f"{PROCESSED_PATH}/OCHEM.Water.Solublity.05.27.2019.smi"
    logging.info(f"Processing {fname}")

    cmpd_list = []
    with open(fname, 'r') as fin:
        for line in fin:
            try:
                pairs = line.rstrip('\n').split(",")
                canon_smiles = canonicalize_smiles(pairs[0])
                logS = float(pairs[1])
                cmpd_list.append((canon_smiles, logS))
            except:
                smiles = pairs[0]
                logging.info(f"Failed to process {smiles}")

    with open(fout, 'w', encoding='ascii') as fo:
        for el in cmpd_list:
            smiles = el[0]
            logS = el[1]
            fo.write(f"{smiles},{logS}\n")

    logging.info(f"Saved {fout}")


def process_A_2019_ADMET_DMPK():
    """ Process A.2019_ADMET_DMPK data.
    Two values are gien in the raw file, so two output files are saved.
    """

    fname = f"{DATA_PATH}/A.2019.ADMET_DMPK.csv"
    fout_1 = f"{PROCESSED_PATH}/A.2019.ADMET_DMPK.SSF.smi"
    fout_2 = f"{PROCESSED_PATH}/A.2019.ADMET_DMPK.CS.smi"
    logging.info(f"Processing {fname}")

    try:
        import pubchempy as pcp
    except ModuleNotFoundError as e:
        print(e)
        return

    with open(fname, 'r') as fin, open(fout_1, 'w') as fout1, open(fout_2, 'w') as fout2:
        fin.readline()
        for line in fin:
            if line.startswith("\""):
                pairs = line.rstrip().split("\"")
                name = pairs[1]
                pairs = pairs[2].split(',')
                logS0_SFF = pairs[0]
                logS0_CS = pairs[2]
            else:
                pairs = line.rstrip().split(',')
                name = pairs[0]
                logS0_SFF = pairs[1]
                logS0_CS = pairs[3]

            name = name.replace('\"', '')
            results = pcp.get_compounds(name, 'name')
            if len(results) > 0:
                isomeric_smiles = results[0].isomeric_smiles
                canon_smiles = canonicalize_smiles(isomeric_smiles)
                fout1.write("{},{}\n".format(canon_smiles, logS0_SFF))
                fout2.write("{},{}\n".format(canon_smiles, logS0_CS))


def process_test_100():
    """
    # SMILES, Interlab.SD, Num.Lit.Sources, Experimental.MP.(*C),log.Poct-water.calc.in.RDKit, log.S0.calc.by.GSE
    """
    fname = f"{TEST_PATH}/set_100.csv"
    fout = f"{TEST_PATH}/test_100.smi"
    fout_gse = f"{TEST_PATH}/test_100.with.gse.smi"
    logging.info(f"Processing {fname}")

    cmpd_list = []
    with open(fname, 'r') as fin:
        for line in fin:
            if not line.startswith("#"):
                pairs = line.rstrip('\n').split(',')
                canon_smiles = canonicalize_smiles(pairs[0])
                logS = float(pairs[5])
                cmpd_list.append((canon_smiles, logS))

    with open(fout, 'w', encoding="ascii") as fo:
        for el in cmpd_list:
            smiles = el[0]
            fo.write(f"{smiles}\n")

    with open(fout_gse, 'w', encoding="ascii") as fo:
        for el in cmpd_list:
            smiles = el[0]
            logS = el[1]
            fo.write(f"{smiles},{logS}\n")

    logging.info(f"Saved {fout}")


def process_test_32():
    """
    # SMILES,
      Interlab.SD,
      Num.Lit.Sources,
      Experimental.MP.(*C),
      log.Poct-water.calc.in.RDKit,
      log.S0.calc.by.GSE
    """
    fname = f"{TEST_PATH}/set_32.csv"
    fout = f"{TEST_PATH}/test_32.smi"
    fout_gse = f"{TEST_PATH}/test_32.with.gse.smi"
    logging.info(f"Processing {fname}")

    cmpd_list = []
    with open(fname, 'r') as fin:
        for line in fin:
            if not line.startswith("#"):
                pairs = line.rstrip('\n').split(',')
                canon_smiles = canonicalize_smiles(pairs[0])
                logS = float(pairs[5])
                cmpd_list.append((canon_smiles, logS))

    with open(fout, 'w', encoding="ascii") as fo:
        for el in cmpd_list:
            smiles = el[0]
            fo.write(f"{smiles}\n")

    with open(fout_gse, 'w', encoding="ascii") as fo:
        for el in cmpd_list:
            smiles = el[0]
            logS = el[1]
            fo.write(f"{smiles},{logS}\n")

    logging.info(f"Saved {fout}")


def unique():
    FILES_LIST = ['BOM.2017.JC.smi',
                  'D.2008.JCIC.smi',
                  'H.2000.JCIC.test1.smi',
                  'H.2000.JCIC.test2.smi',
                  'H.2000.JCIC.train.smi',
                  'HXZ.2004.JCIC.data_set.smi',
                  'HXZ.2004.JCIC.test.smi',
                  'LGG.2008.JCIM.100.smi',
                  'LGG.2008.JCIM.32.smi',
#                  'OCHEM.Water.Solublity.05.27.2019.smi',
                  'POG.2007.JCIM.test.smi',
                  'POG.2007.JCIM.train.smi',
                  'WKH.2007.JCIM.smi',
                  'A.2019.ADMET_DMPK.SSF.smi',
                  'A.2019.ADMET_DMPK.CS.smi'
#                  "WHX.2009.JCIM.Set-001.smi",
#                  "WHX.2009.JCIM.Set-002.smi",
#                  "WHX.2009.JCIM.Set-003.smi",
#                  "WHX.2009.JCIM.Set-004.smi",
#                  "WHX.2009.JCIM.Set-005.smi"
                  ]

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

# Skip the case if the cmpound sign is assigned incorrectly.
            if vmin * vmax < 0:
                continue

# Skip the data if there is too large discrepacy from the
            if (vmax - vmin) > 1.0:
                continue

# Exclude the datapoint if there are some other problems with the measure data
            if sum(np.isinf(values)) != 0:
                continue

# Exclude compounds that are not drug-like
            try:
                if mol_wt(key) > 600.0 or mol_wt(key) < 60.0:
                    continue
            except TypeError as e:
                loggine.error(f"TypeError for {key}: {e}")
                continue

# Save data, comma separated
            fo.write(f"{key};")
            for val in values:
                fo.write(f"{val},")
            fo.write('\n')


def exclude_test():
    '''
    Read all compunds for the trainig set, and exlcude all that are in the test (chellange) set.
    '''

    TEST_32_FILE = f"{TEST_PATH}/test_32.smi"
    TEST_100_FILE = f"{TEST_PATH}/test_100.smi"

    UNIQUE_CMPDS = f"{TRAINING_DIR}/solubility.uniq.smi"

    test_32 = []
    test_100 = []

    with open(TEST_32_FILE, 'r') as fin:
        for line in fin:
            smiles = line.rstrip('\n').split(',')[0]
            test_32.append(smiles)

    with open(TEST_100_FILE, 'r') as fin:
        for line in fin:
            smiles = line.rstrip('\n').split(',')[0]
            test_100.append(smiles)

    unique_dict = {}
    with open(UNIQUE_CMPDS, 'r') as fin:
        for line in fin:
            pairs = line.rstrip('\n').split(';')
            smiles = pairs[0]
            vals_str = pairs[1].split(',')
            vals = [float(v) for v in vals_str[0:-1]]
            unique_dict[smiles] = np.mean(vals)

    TEST_32_FILE_IN_TRAIN = f"{TEST_PATH}/test_32.in-train.smi"
    with open(TEST_32_FILE_IN_TRAIN, 'w') as fo:
        for smiles in test_32:
            if smiles in unique_dict:
                val = unique_dict[smiles]
                fo.write(f"{smiles},{val}\n")

    TEST_100_FILE_IN_TRAIN = f"{TEST_PATH}/test_100.in-train.smi"
    with open(TEST_100_FILE_IN_TRAIN, 'w') as fo:
        for smiles in test_100:
            if smiles in unique_dict:
                val = unique_dict[smiles]
                fo.write(f"{smiles},{val}\n")

    UNIQUE_CMPDS_32 = f"{TRAINING_DIR}/solubility.uniq.no-in-32.smi"
    with open(UNIQUE_CMPDS_32, 'w') as fo:
        for key, value in unique_dict.items():
            if key not in test_32:
                fo.write(f"{key},{value}\n")

    UNIQUE_CMPDS_100 = f"{TRAINING_DIR}/solubility.uniq.no-in-100.smi"
    with open(UNIQUE_CMPDS_100, 'w') as fo:
        for key, value in unique_dict.items():
            if key not in test_100:
                fo.write(f"{key},{value}\n")

    return


def process():
    process_AB_2001_EJPS()
    process_ABB_2000_PR()
    process_BOM_2017_JC()
    process_D_2008_JCIC()
    process_H_2000_JCIC_test1()
    process_H_2000_JCIC_test2()
    process_H_2000_JCIC_train()
    process_HXZ_2004_JCIC_data()
    process_HXZ_2004_JCIC_test()
    process_LGG_2008_JCIM_100()
    process_LGG_2008_JCIM_32()
    process_POG_2007_JCIM_test()
    process_POG_2007_JCIM_train()
    process_WKH_2007_JCIM()
##    process_WHX_2009_JCIM()
    process_A_2019_ADMET_DMPK()
##    process_OCHEM()

    process_test_100()
    process_test_32()


if __name__ == "__main__":
    logging.basicConfig(stream=sys.stdout, level=logging.INFO)

    PROCESSED_DIR = '../data/processed/'
    TRAIN_DIR = '../data/training/'
    TEST_DIR = '../data/test/'

    if not os.path.exists(PROCESSED_DIR):
        os.makedirs(PROCESSED_DIR)

    if not os.path.exists(TRAIN_DIR):
        os.makedirs(TRAIN_DIR)
    
    if not os.path.exists(TEST_DIR):
        os.makedirs(TEST_DIR)

    process()
    unique()
    exclude_test()
