
'''
Confusion matrix of GNPS library searching results
'''
import os
import time

import numpy as np
import pandas as pd
from collections import Counter
from tqdm import tqdm, trange
from rdkit import Chem
from rdkit.Chem import rdFMCS

def MSC(mol1, mol2):

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

if __name__ == '__main__':
    t = time.time()

    # Loading GNPS library searching result
    gnps_file = './data/LibSearch_top100_0.1.csv'
    gnps_df = pd.read_csv(gnps_file,encoding = 'gbk')

    # Categorize into two classes: above thresholds (at) 和 under thresholds (ut)
    scans = list(Counter(gnps_df['#Scan#'])) # Get a list of feature ids
    threshold = 0.7
    mps = 5
    df_at = gnps_df[(gnps_df['MQScore'] >= threshold)&(gnps_df['SharedPeaks'] >= mps)].reset_index(drop=True)
    scans_at = list(Counter(df_at['#Scan#']))

    df_ut = gnps_df[(gnps_df['MQScore'] < threshold)|(gnps_df['SharedPeaks'] < mps)].reset_index(drop=True)
    scans_ut = [x for x in list(Counter(df_ut['#Scan#'])) if x not in scans_at]

    index_to_keep_at = []
    index_to_keep_ut = []

    for i in tqdm(scans_at,total = len(scans_at)):
        slice_df_at = df_at[df_at['#Scan#'] == i]
        max_index_at = slice_df_at['MQScore'].idxmax()
        index_to_keep_at.append(max_index_at)

    for i in tqdm(scans_ut, total=len(scans_ut)):
        slice_df_ut = df_ut[df_ut['#Scan#'] == i]
        max_index_ut = slice_df_ut['MQScore'].idxmax()
        index_to_keep_ut.append(max_index_ut)

    df_at = df_at.loc[index_to_keep_at].reset_index(drop=True)
    df_ut = df_ut.loc[index_to_keep_ut].reset_index(drop=True)

    ## Check if Smiles in edb is missing'''
    # df_at.to_csv('demo_at.csv',index=None)
    # df_ut.to_csv('demo_ut.csv',index=None)

    # MCS calculating
    df_at['mcs'] = np.nan
    for i in trange(len(df_at)):
        smile1 = df_at.loc[i,'Smiles']
        mol1 = Chem.MolFromSmiles(smile1)
        smile2 = df_at.loc[i,'smiles']
        mol2 = Chem.MolFromSmiles(smile2)
        df_at.loc[i,'mcs'] = MSC(mol1,mol2)
    # df_at.to_csv('GNPS_LibSearch_at.csv') # Save the temp dataframe

    df_ut['mcs'] = np.nan
    for i in trange(len(df_ut)):
        smile1 = df_ut.loc[i,'Smiles']
        mol1 = Chem.MolFromSmiles(smile1)
        smile2 = df_ut.loc[i,'smiles']
        mol2 = Chem.MolFromSmiles(smile2)
        try:
            df_ut.loc[i, 'mcs'] = MSC(mol1, mol2)
        except:
            print(i,smile1,smile2)
    # df_ut.to_csv('GNPS_LibSearch_ut.csv') # Save the temp dataframe

    # df_at = pd.read_csv('./data/GNPS_LibSearch_at.csv')
    # df_ut = pd.read_csv('./data/GNPS_LibSearch_ut.csv')

    # Data for visualization (Figure 2E)
    scans_at = list(Counter(df_at['#Scan#']))
    scans_ut = list(Counter(df_ut['#Scan#']))
    max_at = []
    for i in tqdm(scans_at, total=len(scans_at)):
        df_slice = df_at[df_at['#Scan#'] == i]
        max_at.append(max(df_slice['mcs']))
    print('FP:', sum(i < 0.35 for i in max_at)
          , 'TP', sum(i >= 0.7 for i in max_at))

    max_ut = []
    for i in tqdm(scans_ut, total=len(scans_ut)):
        df_slice = df_ut[df_ut['#Scan#'] == i]
        max_ut.append(max(df_slice['mcs']))
    print('TN:', sum(i < 0.35 for i in max_ut)
          , 'FN', sum(i >= 0.7 for i in max_ut))

    print(f'Finished in {(time.time() - t) / 60:.2f} min')
