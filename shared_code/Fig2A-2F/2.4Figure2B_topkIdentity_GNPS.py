# -*- coding: utf-8 -*-
# @Time :2023/6/18 13:57
# @Auther :Yuwenchao
# @Software : PyCharm
'''
计算GNPS比对结果与已知化合物的tanimoto和MCS
然后输出 topK identity 的数量和详细情况

'''
import os
import time

import pandas as pd
import numpy as np

from rdkit import Chem
from rdkit.Chem import rdFMCS
from rdkit.Chem import AllChem
from rdkit.Chem import DataStructs

from tqdm import trange
from collections import Counter

def MCS(mol1, mol2):
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
    os.chdir('/Users/hehe/Desktop/HiTS/MNA/#Fig2A-2F/GNPS_LibSearch_no_analog')

    '''data importing'''
    file = 'LibSearch_top100_0.1.tsv' # GNPS matching results
    df = pd.read_csv(file,sep='\t',encoding='gbk')
    print(df.columns)
    df['tanimoto'] = np.nan
    df['mcs'] = np.nan

    '''计算ecfp4_1024bit'''
    for i in trange(len(df)):
        gnps_smile = df.loc[i,'Smiles']
        gnps_mol = Chem.MolFromSmiles(gnps_smile)
        smile = df.loc[i, 'smiles']
        mol = Chem.MolFromSmiles(smile)
        try:
            gnps_fp = AllChem.GetMorganFingerprintAsBitVect(gnps_mol, radius = 2, nBits=1024)
            fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius = 2, nBits=1024)
            sim = DataStructs.TanimotoSimilarity(gnps_fp,fp)
            mcs = MCS(mol,gnps_mol)
        except:print(i)
        df.loc[i,'tanimoto'] = sim
        df.loc[i, 'mcs'] = mcs


    # '''Get identity in topK'''
    # top1, top2, top3, top4, top5, top10, top20 =0, 0, 0, 0, 0, 0, 0
    # scans = list(Counter(df['#Scan#']))
    # index_edb_to_keep1,index_edb_to_keep5,index_edb_to_keep20 = [],[],[]
    # for scan in scans:
    #     df_slice = df[df['#Scan#']==scan]
    #     df_slice_sorted = df_slice.sort_values('MQScore', ascending=False)
    #     '''统计topK的identity数量'''
        # top_1 = df_slice_sorted.head(1)['tanimoto'].values.tolist()
        # top_5 = df_slice_sorted.head(5)['tanimoto'].values.tolist()
        # top_20 = df_slice_sorted.head(20)['tanimoto'].values.tolist()
        # if 1.0 in top_1:
        #     top1+=1
        # if 1.0 in top_5:
        #     top5+=1
        # if 1.0 in top_20:
        #     top20+=1
        # '''输出'''
    #     index_top1 = df_slice_sorted.head(1).index
    #     index_edb_to_keep1.extend(index_top1)
    #
    #     index_top5 = df_slice_sorted.head(5).index
    #     index_edb_to_keep5.extend(index_top5)
    #
    #     index_top20 = df_slice_sorted.head(20).index
    #     index_edb_to_keep20.extend(index_top20)
    #
    # print(top1, top5, top20)
    #
    # index_edb_to_keep1 = list(set(index_edb_to_keep1))
    # index_edb_to_keep5 = list(set(index_edb_to_keep5))
    # index_edb_to_keep20 = list(set(index_edb_to_keep20))
    # edb_df_top1 = df.loc[index_edb_to_keep1]
    # edb_df_top5 = df.loc[index_edb_to_keep5]
    # edb_df_top20 = df.loc[index_edb_to_keep20]
    # edb_df_top1.to_csv('gnps_top1.csv')
    # edb_df_top5.to_csv('gnps_top5.csv')
    # edb_df_top20.to_csv('gnps_top20.csv')


    print(f'Finished in {(time.time() - t) / 60:.2f} min')
