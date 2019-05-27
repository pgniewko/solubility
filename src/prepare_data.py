import sys
import logging

import numpy as np
from utils import canonicalize_smiles, mol_wt
from rdkit import Chem

DATA_PATH="/Users/pawel/Projects/solubility/data/raw"
PROCESSED_PATH="/Users/pawel/Projects/solubility/data/processed"
TEST_PATH="/Users/pawel/Projects/solubility/data/test"

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
                cmpd_list.append((canon_smiles,logS))

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
            cmpd_list.append((canon_smiles,logS))

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
                cmpd_list.append((canon_smiles,logS))

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
                pairs = line.rstrip('\n').split()
                canon_smiles = canonicalize_smiles(pairs[6])
                logS = float(pairs[3])
                cmpd_list.append((canon_smiles,logS))

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
                pairs = line.rstrip('\n').split()
                canon_smiles = canonicalize_smiles(pairs[6])
                logS = float(pairs[3])
                cmpd_list.append((canon_smiles,logS))
                
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
                pairs = line.rstrip('\n').split()
                canon_smiles = canonicalize_smiles(pairs[6])
                logS = float(pairs[3])
                cmpd_list.append((canon_smiles,logS))

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
        cnt=0
        for line in fin:
            if not line.startswith("#"):
                pairs = line.rstrip('\n').split()
                canon_smiles = canonicalize_smiles(pairs[0])
                logS = float(pairs[2])
                cmpd_list.append((canon_smiles,logS))
                cnt+=1

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
        cnt=0
        for line in fin:
            if not line.startswith("#"):
                pairs = line.rstrip('\n').split(',')
                canon_smiles = canonicalize_smiles(pairs[0])
                logS = float(pairs[2])
                cmpd_list.append((canon_smiles,logS))
                cnt+=1

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
        cnt=0
        for line in fin:
            if not line.startswith("#"):
                pairs = line.rstrip('\n').split(',')
                canon_smiles = canonicalize_smiles(pairs[0])
                mg = float(pairs[1]) / 1000
                mw = mol_wt(canon_smiles)
                logS = np.log10(mg/mw)
                cmpd_list.append((canon_smiles,logS))
                cnt+=1

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
        cnt=0
        for line in fin:
            if not line.startswith("#"):
                pairs = line.rstrip('\n').split(";")
                canon_smiles = canonicalize_smiles(pairs[1].strip("\""))
                logS = float(pairs[4])
                cmpd_list.append((canon_smiles,logS))
                cnt+=1

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
        cnt=0
        for line in fin:
            if not line.startswith("#"):
                pairs = line.rstrip('\n').split(";")
                canon_smiles = canonicalize_smiles(pairs[1].strip("\""))
                logS = float(pairs[4])
                cmpd_list.append((canon_smiles,logS))
                cnt+=1

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
        cnt=0
        for line in fin:
            if not line.startswith("#"):
                pairs = line.rstrip('\n').split(',')
                canon_smiles = canonicalize_smiles(pairs[0])
                logS = float(pairs[5])
                cmpd_list.append((canon_smiles ,logS))
                cnt+=1

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
    # SMILES, Interlab.SD, Num.Lit.Sources, Experimental.MP.(*C),log.Poct-water.calc.in.RDKit, log.S0.calc.by.GSE
    """
    fname = f"{TEST_PATH}/set_32.csv"
    fout = f"{TEST_PATH}/test_32.smi"
    fout_gse = f"{TEST_PATH}/test_32.with.gse.smi"
    logging.info(f"Processing {fname}")

    cmpd_list = []
    with open(fname, 'r') as fin:
        cnt=0
        for line in fin:
            if not line.startswith("#"):
                pairs = line.rstrip('\n').split(',')
                canon_smiles = canonicalize_smiles(pairs[0])
                logS = float(pairs[5])
                cmpd_list.append((canon_smiles ,logS))
                cnt+=1

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


def main():
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

    process_test_100()
    process_test_32()

if __name__ == "__main__":
    logging.basicConfig(stream=sys.stdout, level=logging.DEBUG)
    main()
