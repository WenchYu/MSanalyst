
'''
Calculating MCS similarity
Calculating theoretical mass of adducts
    smiles -> formula -> exact mass of adducts
'''
import re
import pandas as pd
import numpy as np
from tqdm import trange
from rdkit import Chem
from rdkit.Chem import rdFMCS
from rdkit.Chem import rdMolDescriptors


def MCS(mol1, mol2):
    '''
    Calculating maximum common substructure (MCS) similarity
    from [literature](https://www.nature.com/articles/s41467-022-30118-9)
    '''
    mcs = rdFMCS.FindMCS([mol1, mol2]
                         , bondCompare=rdFMCS.BondCompare.CompareOrder
                         , atomCompare=rdFMCS.AtomCompare.CompareAny
                         , maximizeBonds = False
                         , ringMatchesRingOnly=False
                         , matchValences=False
                         , timeout=10
                         )
    mcs_num_bonds = mcs.numBonds
    mol1_num_bonds = mol1.GetNumBonds()
    mol2_num_bonds = mol2.GetNumBonds()
    similarity = mcs_num_bonds / ((mol1_num_bonds + mol2_num_bonds) - mcs_num_bonds)
    return similarity

def Smile2Formula(smile):
    '''
    Converting smile to chemical formula
    '''
    mol = Chem.MolFromSmiles(smile)
    if mol is None:
        return "Invalid SMILES"
    mol_with_h = Chem.AddHs(mol)
    formula = rdMolDescriptors.CalcMolFormula(mol_with_h)
    return formula

class MyChemInfo():
    @staticmethod
    def AtomicWeight(element: str) -> float:
        """
        Monoisotopic mass of each element
        """
        if len(element) > 2:
            return None
        return {
            "H": 1.007825,
            "C": 12.0000,
            "N": 14.003074,
            "O": 15.994915,
            "F": 18.000938,
            "He": 4.002602,
            "Li": 6.94,
            "Be": 9.0121831,
            "B": 10.012937,
            "Cl": 34.968853,
            "Br": 78.918337,
            "Ne": 20.1797,
            "Na": 22.98976928,
            "Mg": 24.305,
            "Al": 26.9815385,
            "Si": 27.976927,
            "P": 30.973770,
            "S": 31.972071,
            "I": 126.904473,
            "Ar": 39.948,
            "K": 39.0983,
            "Ca": 40.078,
            "Sc": 44.955908,
            "Ti": 47.867,
            "V": 50.9415,
            "Cr": 51.9961,
            "Mn": 54.938044,
            "Fe": 55.845,
            "Co": 58.933194,
            "Ni": 58.6934,
            "Cu": 63.546,
            "Zn": 65.38,
            "Ga": 69.723,
            "Ge": 72.63,
            "As": 74.921595,
            "Se": 73.922476,
            "Kr": 83.798,
            "Rb": 85.4678,
            "Sr": 87.62,
            "Y": 88.90584,
            "Zr": 91.224,
            "Nb": 92.90637,
            "Mo": 95.95,
            "Ru": 101.07,
            "Rh": 102.9055,
            "Pd": 106.42,
            "Ag": 107.8682,
            "Cd": 112.414,
            "In": 114.818,
            "Sn": 118.71,
            "Sb": 121.76,
            "Te": 127.6,
            "Xe": 131.293,
            "Cs": 132.90545196,
            "Ba": 137.327,
            "La": 138.90547,
            "Ce": 140.116,
            "Pr": 140.90766,
            "Nd": 144.242,
            "Sm": 150.36,
            "Eu": 151.964,
            "Gd": 157.25,
            "Tb": 158.92535,
            "Dy": 162.5,
            "Ho": 164.93033,
            "Er": 167.259,
            "Tm": 168.93422,
            "Yb": 173.054,
            "Lu": 174.9668,
            "Hf": 178.49,
            "Ta": 180.94788,
            "W": 183.84,
            "Re": 186.207,
            "Os": 190.23,
            "Ir": 192.217,
            "Pt": 195.084,
            "Au": 196.966569,
            "Hg": 200.592,
            "Tl": 204.38,
            "Pb": 207.2,
            "Bi": 208.9804,
            "Th": 232.0377,
            "Pa": 231.03588,
            "U": 238.02891,
            "Tc": 0,
            "Pm": 0,
            "Po": 0,
            "At": 0,
            "Rn": 0,
            "Fr": 0,
            "Ra": 0,
            "Ac": 0,
            "Np": 0,
            "Pu": 0,
            "Am": 0,
            "Cm": 0,
            "Bk": 0,
            "Cf": 0,
            "Es": 0,
            "Fm": 0,
            "Md": 0,
            "No": 0,
            "Lr": 0,
            "Rf": 0,
            "Db": 0,
            "Sg": 0,
            "Bh": 0,
            "Hs": 0,
            "Mt": 0,
            "Ds": 0,
            "Rg": 0,
            "Cn": 0,
            "Fl": 0,
            "Lv": 0}.get(element, 0.000)

    @staticmethod
    def MolWt(formula: str) -> float:
        '''
        Calculating the exact mass from a chemical formula
        :param formula: Chemical formular in str format
        '''
        regStr = "([A-Z]{1}[a-z]{0,1})([0-9]{0,3})"
        MatchList = re.findall(regStr, formula)
        cntMatchList = len(MatchList)
        i = 0
        mW = 0.000
        while i < cntMatchList:
            eleName = MatchList[i][0]
            eleCount = int(MatchList[i][1]) if len(MatchList[i][1]) > 0 else 1
            aw = MyChemInfo.AtomicWeight(eleName)
            if (aw == 0):
                return 0
            mW += aw * eleCount
            i = i + 1
        return mW

    @staticmethod
    def Adduct(adduct:str)-> float:
        regStr = "([A-Z]{1}[a-z]{0,1})([0-9]{0,3})"
        MatchList = re.findall(regStr, adduct)
        for i in range(len(MatchList)):
            if MatchList[i][0] != 'M':
                return MyChemInfo.AtomicWeight(MatchList[i][0])

if __name__ == '__main__':
    '''Loading dataset'''
    MS1_library_file = '../msdb/isdbMS1.csv' # Example csv file, containing columns['smiles']
    MS1_library_df = pd.read_csv(MS1_library_file,index_col=None)

    '''Calculating the theoretical values of all adducts'''
    MS1_library_df['formula'] = np.nan
    MS1_library_df['exactmass'] = np.nan
    for i in trange(len(MS1_library_df.index)):
        try:
            smile = MS1_library_df.smiles[i]
            MS1_library_df.loc[i,'formula'] = Smile2Formula(smile)
            MS1_library_df.loc[i,'exactmass'] = MyChemInfo.MolWt(MS1_library_df.formula[i])
            MS1_library_df.loc[i,'m+h'] = MS1_library_df.loc[i,'exactmass'] + 1.007276
            MS1_library_df.loc[i,'m+nh4'] = MS1_library_df.loc[i,'exactmass'] + 18.033823
            MS1_library_df.loc[i,'m+na'] = MS1_library_df.loc[i,'exactmass'] + 22.989218
        except:
            pass










