#!/usr/bin/env python

import sys
from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen, Lipinski
from rdkit.Chem import PandasTools
import pandas as pd
from sklearn.linear_model import LinearRegression
from collections import namedtuple


class ESOLCalculator:
    def __init__(self):
        self.aromatic_query = Chem.MolFromSmarts("a")
        self.Descriptor = namedtuple("Descriptor", "mw logp rotors ap")
        self._coef = {"MW":0.0, "LogP":0.0, "RB":0.0, "AP":0.0}
        self._intercept = 0.0


    def fit(self, smiles_list, logS_list):
        descriptor_cols = ["MW", "Logp", "Rotors", "AP"]
        X = []
        y = []
        for i, smiles in enumerate(smiles_list):
            mol = Chem.MolFromSmiles(smiles)
            (mw, logp, rotors, ap) = self._calc_esol_descriptors(mol)
            X.append( [mw, logp, rotors, ap] )
            y.append(logS_list[i])

        model = LinearRegression()
        model.fit(X, y)
        print("score:", model.score(X, y))
        coefficient_dict = {}
        #self._intercept    =  0.2612 #model.intercept_
        #self._coef["MW"]   = -0.0066 #model.coef_[0]
        #self._coef["LogP"] = -0.7417 #model.coef_[1]
        #self._coef["RB"]   =  0.0035 #model.coef_[2]
        #self._coef["AP"]   = -0.4262 #model.coef_[3]
        self._intercept    = model.intercept_
        self._coef["MW"]   = model.coef_[0]
        self._coef["LogP"] = model.coef_[1]
        self._coef["RB"]   = model.coef_[2]
        self._coef["AP"]   = model.coef_[3]
        print("_score:", self.score(smiles_list, y))
        
        #mw': -0.0066138847738667125, 'logp': -0.7416739523408995, 'rotors': 0.003451545565957996, 'ap': -0.42624840441316975

    def predict(self, smiles_list):
        ypred = []
        for smiles in smiles_list:
            mol = Chem.MolFromSmiles(smiles)
            (mw, logp, rotors, ap) = self._calc_esol_descriptors(mol)
            y_val  = self._intercept
            y_val += self._coef["MW"] * mw
            y_val += self._coef["LogP"] * logp
            y_val += self._coef["RB"] * rotors
            y_val += self._coef["AP"] * ap 
            ypred.append(y_val)

        return ypred


    def score(self, smiles_list, y_measured):
        ypred = self.predict(smiles_list)
        tot = len(y_measured)

        mse = 0.0
        mae = 0.0
        for i, yval in enumerate(ypred):
            mse += (yval - y_measured[i])**2.0
            mae += abs(yval - y_measured[i])
        
        mse /= tot
        mae /= tot

        return mse, mae


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
        :return: named tuple with descriptor values
        """
        mw = Descriptors.MolWt(mol)
        logp = Crippen.MolLogP(mol)
        rotors = Lipinski.NumRotatableBonds(mol)
        ap = self._calc_ap(mol)
        return (mw, logp, rotors, ap)

        #return self.Descriptor(mw=mw, logp=logp, rotors=rotors, ap=ap)

#    def calc_esol_orig(self, mol):
#        """
#        Original parameters from the Delaney paper, just here for comparison
#        :param mol: input molecule
#        :return: predicted solubility
#        """
#        # just here as a reference don't use this!
#        intercept = 0.16
#        coef = {"logp": -0.63, "mw": -0.0062, "rotors": 0.066, "ap": -0.74}
#        desc = self.calc_esol_descriptors(mol)
#        esol = intercept + coef["logp"] * desc.logp + coef["mw"] * desc.mw + coef["rotors"] * desc.rotors \
#               + coef["ap"] * desc.ap
#        return esol


#    def calc_esol(self, mol):
#        """
#        Calculate ESOL based on descriptors in the Delaney paper, coefficients refit for the RDKit using the
#        routine refit_esol below
#        :param mol: input molecule
#        :return: predicted solubility
#        """
#        intercept = 0.26121066137801696
#        coef = {'mw': -0.0066138847738667125, 'logp': -0.7416739523408995, 'rotors': 0.003451545565957996, 'ap': -0.42624840441316975}
#        desc = self.calc_esol_descriptors(mol)
#        esol = intercept + coef["logp"] * desc.logp + coef["mw"] * desc.mw + coef["rotors"] * desc.rotors \
#               + coef["ap"] * desc.ap
#        return esol


#def test_on_dls_100():
#    """
#    Test on the dls_100 dataset from the University of St Andrews
#    https://doi.org/10.17630/3a3a5abc-8458-4924-8e6c-b804347605e8
#    This function downloads the Excel spreadhsheet from St Andrews, does the
#    comparison.  Note that this does not remove the structures in DLS-100 that are
#    also in the training set.
#    :return: None
#    """
#    esol_calculator = ESOLCalculator()
#    df = pd.read_excel("https://risweb.st-andrews.ac.uk/portal/files/251452157/dls_100.xlsx")
#    df.dropna(inplace=True)
#    PandasTools.AddMoleculeColumnToFrame(df, 'SMILES', 'Molecule', includeFingerprints=False)
#    res = []
#    for mol, name, logS in df[["Molecule", "Chemical name", "LogS exp (mol/L)"]].values:
#        res.append([name, logS, esol_calculator.calc_esol(mol)])
#    df = pd.DataFrame(res, columns=["Name", "Experimental LogS", "ESol LogS"])
#    df.to_csv("dls_100.csv", index=False)


#def add_esol_descriptors_to_dataframe(df,smiles_col="SMILES",name_col="Compound ID"):
#    """
#    Add ESOL descriptors to a Pandas dataframe
#    :param df: input dataframe
#    :param smiles_col: column in the dataframe with SMILES
#    :param name_col: column with the molecule names
#    :return: dataframe and list of columns that were added
#    """
#    esol_calculator = ESOLCalculator()
#    PandasTools.AddMoleculeColumnToFrame(df, smiles_col, 'Molecule', includeFingerprints=False)
#    result_list = []
#    for name, mol in df[[name_col, "Molecule"]].values:
#        desc = esol_calculator.calc_esol_descriptors(mol)
#        result_list.append([name,desc.mw,desc.logp,desc.rotors,desc.ap])
#    result_df = pd.DataFrame(result_list)
#    descriptor_cols = ["MW", "Logp", "Rotors", "AP"]
#    result_df.columns = [name_col] + descriptor_cols
#    df = df.merge(result_df, on=name_col)
#    return df, descriptor_cols


#def refit_esol(input_file_name, truth_col):
#    """
#    Refit the parameters for ESOL using multiple linear regression
#    Prints parameters that can be pasted into the calc_esol function
#    :input_file_name: input file
#    :truth_col: column with the experimental value
#    :return: None
#    """
#    df = pd.read_csv(input_file_name)
#    df, descriptor_cols = add_esol_descriptors_to_dataframe(df)
#    x = df[descriptor_cols]
#    y = df[[truth_col]]
#
#    model = LinearRegression()
#    model.fit(x, y)
#    coefficient_dict = {}
#    for name, coef in zip(descriptor_cols, model.coef_[0]):
#        coefficient_dict[name.lower()] = coef
#    print("intercept = ", model.intercept_[0])
#    print("coef =", coefficient_dict)


#def demo(truth_col):
#    """
#    Read the csv file from the Delaney supporting information, calculate ESOL using the data file from the
#    supporting material using the original and refit coefficients.  Write a csv file comparing with experiment.
#    :return: None
#    """
#    esol_calculator = ESOLCalculator()
#    df = pd.read_csv("delaney.csv")
#    PandasTools.AddMoleculeColumnToFrame(df, 'SMILES', 'Molecule', includeFingerprints=False)
#    res = []
#    for mol, val in df[["Molecule", truth_col]].values:
#        res.append([val, esol_calculator.calc_esol(mol), esol_calculator.calc_esol_orig(mol)])
#    output_df = pd.DataFrame(res, columns=["Experiment", "ESOL Current", "ESOL Original"])
#    output_df.to_csv('validation.csv',index=False)
#

#def main(smiles_list, logS_list):
#
#
#    truth_col = "measured log(solubility:mol/L)"
#    refit_esol("delaney.csv", truth_col)
#    demo(truth_col)


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


if __name__ == "__main__":

    smiles_list, logS_list = get_training_data(sys.argv[1])
    esol_calculator = ESOLCalculator()
    esol_calculator.fit(smiles_list, logS_list)
    y_pred = esol_calculator.predict(smiles_list)
    esol_calculator.score(smiles_list, logS_list)

    import matplotlib.pyplot as plt
    plt.plot(logS_list, y_pred, 'o', alpha=0.01)
    plt.plot([-10,2],[-10,2], '--')
    plt.show()

