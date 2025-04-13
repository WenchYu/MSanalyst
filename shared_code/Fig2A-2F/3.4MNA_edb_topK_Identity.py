# -*- coding: utf-8 -*-
# @Time :2023/6/18 16:42
# @Auther :Yuwenchao
# @Software : PyCharm
'''
使用 tanimoto 来判断 identity 的 molecule
'''
import os
import time

import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import DataStructs

from tqdm import trange
from collections import Counter

if __name__ == '__main__':
    t = time.time()
    os.chdir('/Users/hehe/Desktop/HTS/MNA/#std_mix/match_results_mna_10ppm')

    '''Data importing'''
    edb_file = 'E_MS1match_std_quant.csv'
    isdb_file = 'IS_MS1match_std_quant.csv'

    edb_df = pd.read_csv(edb_file)
    isdb_df = pd.read_csv(isdb_file)
    print(edb_df.columns)

    '''计算1024_ecfp4'''
    for i in trange(len(edb_df)):
        try:
            edb_smile = edb_df.loc[i, 'match_smiles']
            edb_mol = Chem.MolFromSmiles(edb_smile)
            smile = edb_df.loc[i, 'smiles']
            mol = Chem.MolFromSmiles(smile)
            gnps_fp = AllChem.GetMorganFingerprintAsBitVect(edb_mol, radius=2, nBits=1024)
            fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=1024)
            sim = DataStructs.TanimotoSimilarity(gnps_fp, fp)
            edb_df.loc[i, 'tanimoto'] = sim
        except:
            pass

    '''<Check>'''
    edb_df.to_csv('test.csv')

    '''Get identity in topK'''
    top1, top2, top3, top4, top5, top10, top20 = 0, 0, 0, 0, 0, 0, 0
    edb_scans = list(Counter(edb_df['row ID']))
    index_edb_to_keep1, index_edb_to_keep5, index_edb_to_keep20 = [], [], []
    for edb_scan in edb_scans:
        edb_df_slice = edb_df[edb_df['row ID'] == edb_scan]
        edb_df_slice_sorted = edb_df_slice.sort_values('pp',ascending = False)
        '''统计topK的identity数量'''
        # top_1 = edb_df_slice_sorted.head(1)['tanimoto'].values.tolist()
        # top_5 = edb_df_slice_sorted.head(5)['tanimoto'].values.tolist()
        # top_20 = edb_df_slice_sorted.head(20)['tanimoto'].values.tolist()
        # if 1.0 in top_1:
        #     top1 += 1
        # if 1.0 in top_5:
        #     top5 += 1
        # if 1.0 in top_20:
        #     top20 += 1

        index_top1 = edb_df_slice_sorted.head(1).index
        index_edb_to_keep1.extend(index_top1)

        index_top5 = edb_df_slice_sorted.head(5).index
        index_edb_to_keep5.extend(index_top5)

        index_top20 = edb_df_slice_sorted.head(20).index
        index_edb_to_keep20.extend(index_top20)

    print(top1, top5, top20)

    index_edb_to_keep1 = list(set(index_edb_to_keep1))
    index_edb_to_keep5 = list(set(index_edb_to_keep5))
    index_edb_to_keep20 = list(set(index_edb_to_keep20))
    edb_df_top1 = edb_df.loc[index_edb_to_keep1]
    edb_df_top5 = edb_df.loc[index_edb_to_keep5]
    edb_df_top20 = edb_df.loc[index_edb_to_keep20]
    edb_df_top1.to_csv('edb_pp_top1.csv')
    edb_df_top5.to_csv('edb_pp_top5.csv')
    edb_df_top20.to_csv('edb_pp_top20.csv')





    print(f'Finished in {(time.time() - t) / 60:.2f} min')
