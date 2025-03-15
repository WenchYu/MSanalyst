# -*- coding: utf-8 -*-
# @Time :2023/6/10 10:24
# @Auther :Yuwenchao
# @Software : PyCharm
'''
Confusion matrix of msanalyst results
'''
import os
import time

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import rdFMCS
from collections import Counter
from tqdm import tqdm,trange

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
    '''Importing data'''
    isdb_file = f'data/npMS1match_std_quant.csv'
    edb_file = f'data/edbMS1match_std_quant.csv'
    isdb_df = pd.read_csv(isdb_file)
    edb_df = pd.read_csv(edb_file)

    '''Take the maximum value of pp and pair_similarity for each row.'''
    isdb_df['pair_similarity'] = np.nan # set the colums of isdb_match_file the same as edb
    isdb_df['mps'] = np.nan
    isdb_df['pp'] = np.nan
    for j in trange(len(isdb_df)):
        pairs = [
                (isdb_df.loc[j, 'mps0'], isdb_df.loc[j, 'pair_similarity0'], isdb_df.loc[j, 'pp0']),
                (isdb_df.loc[j, 'mps1'], isdb_df.loc[j, 'pair_similarity1'], isdb_df.loc[j, 'pp1']),
                (isdb_df.loc[j, 'mps2'], isdb_df.loc[j, 'pair_similarity2'], isdb_df.loc[j, 'pp2'])]

        '''Find the (mps, pairsimilarity) pair with the maximum pairsimilarity value'''
        max_sim_pair = max(pairs, key=lambda x: x[1]) # [0] mps [1] pair_similarity [2] pp
        # Retrieve the maximum mps and pair_similarity values
        isdb_df.loc[j, 'mps'] = max_sim_pair[0]
        isdb_df.loc[j, 'pair_similarity'] = max_sim_pair[1]

        '''Find the (mps, pp) pair with the maximum pp value'''
        max_pp_pair = max(pairs, key=lambda x: x[2])  # [0] mps [1] pair_similarity [2] pp
        # Retrieve the maximum mps and pair_similarity values
        isdb_df.loc[j, 'pp'] = max_sim_pair[2]

    columns_to_merge = ['row ID', 'smiles','match_smiles', 'pair_similarity', 'mps', 'pp']
    merged_data = pd.concat(
        [isdb_df[columns_to_merge], edb_df[columns_to_merge]], axis=0)

    '''Taking the maximum values of pp and pair similarity while retaining the corresponding indices, and then calculate their MCS'''
    scans = list(Counter(isdb_df['row ID']))
    edb_index_to_keep = []
    isdb_index_to_keep = []
    for i in tqdm(scans, total = len(scans)):
        slice_edb_df = edb_df[edb_df['row ID'] == i]
        max_sim_index1 = slice_edb_df['pair_similarity'].idxmax()
        max_pp_index1 = slice_edb_df['pp'].idxmax()
        edb_index_to_keep.append(max_sim_index1)
        edb_index_to_keep.append(max_pp_index1)

        slice_isdb_df = isdb_df[isdb_df['row ID'] == i]
        max_sim_index2 = slice_isdb_df['pair_similarity'].idxmax()
        max_pp_index2 = slice_isdb_df['pp'].idxmax()
        isdb_index_to_keep.append(max_sim_index2)
        isdb_index_to_keep.append(max_pp_index2)

    edb_index_to_keep = list(set(edb_index_to_keep)) # 对index去重
    isdb_index_to_keep = list(set(isdb_index_to_keep))
    edb_filtered = edb_df.loc[edb_index_to_keep].reset_index(drop=True)
    isdb_filtered = isdb_df.loc[isdb_index_to_keep].reset_index(drop=True)

    '''计算MCS相似度'''
    edb_filtered['mcs'] = np.nan
    isdb_filtered['mcs'] = np.nan
    for i in trange(len(edb_filtered)):
        smile1 = edb_filtered.loc[i,'smiles']
        mol1 = Chem.MolFromSmiles(smile1)
        smile2 = edb_filtered.loc[i,'match_smiles']
        mol2 = Chem.MolFromSmiles(smile2)
        edb_filtered.loc[i,'mcs'] = MSC(mol1,mol2)

    for j in trange(len(isdb_filtered)):
        smile1 = isdb_filtered.loc[j,'smiles']
        mol1= Chem.MolFromSmiles(smile1)
        smile2 = isdb_filtered.loc[j,'match_smiles']
        mol2 = Chem.MolFromSmiles(smile2)
        isdb_filtered.loc[j,'mcs'] =MSC(mol1,mol2)

    # edb_filtered.to_csv('edb_fil.csv',index=None) # MCS calculation takes a long time to run
    # isdb_filtered.to_csv('isdb_fil.csv',index=None)
    # edb_filtered = pd.read_csv('edb_fil.csv')
    # isdb_filtered = pd.read_csv('isdb_fil.csv')

    '''
    Merge the results from edb and isdb
    Categrize them into 2 class(selected or unselected) based on the thresholds
    Counts based on the MCS values
    '''

    columns_to_merge = ['row ID', 'pair_similarity', 'mps', 'pp','mcs']
    merged_data = pd.concat(
        [edb_filtered[columns_to_merge], isdb_filtered[columns_to_merge]], axis=0)

    scans = list(Counter(edb_filtered['row ID']))
    threshold = 0.7
    mps = 5
    df_selected = merged_data[((merged_data['pair_similarity'] >= threshold) & (merged_data['mps']>=mps))|(
            (merged_data['pp']>=threshold)&(merged_data['mps']>=mps))]
    df_unselected = merged_data[~(((merged_data['pair_similarity'] >= threshold) & (merged_data['mps'] >= mps)) | (
                (merged_data['pp'] >= threshold) & (merged_data['mps'] >= mps)))]

    scans_selected = list(Counter(df_selected['row ID']))
    max_selected= []
    for i in tqdm(scans_selected,total = len(scans_selected)):
        df_slice = df_selected[df_selected['row ID'] == i]
        max_selected.append(max(df_slice['mcs']))
    print('FP:',sum(i < 0.35 for i in max_selected)
          ,'TP',sum(i >= 0.7 for i in max_selected))

    scans_unselected = [x for x in list(Counter(df_unselected['row ID'])) if x not in scans_selected]
    max_unselected = []
    for i in tqdm(scans_unselected, total=len(scans_unselected)):
        df_slice = df_unselected[df_unselected['row ID'] == i]
        max_unselected.append(max(df_slice['mcs']))
    print('TN:', sum(i < 0.35 for i in max_unselected)
          , 'FN', sum(i >= 0.7 for i in max_unselected))

    print(f'Finished in {(time.time() - t) / 60:.2f} min')
