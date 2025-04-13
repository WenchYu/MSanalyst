# -*- coding: utf-8 -*-
# @Time :2023/6/18 16:42
# @Auther :Yuwenchao
# @Software : PyCharm
'''

'''
import os
import time

import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import DataStructs

from tqdm import trange
from collections import Counter

if __name__ == '__main__':
    t = time.time()
    os.chdir('/Users/hehe/Desktop/HTS/MNA/#std_mix/match_results_mna_10ppm')

    '''Data importing'''
    isdb_file = 'IS_MS1match_std_quant.csv'

    isdb_df = pd.read_csv(isdb_file)
    print(isdb_df.columns)
    '''处理一下ISDB比对结果，取每一行pp和pair_simialrity的最大值'''
    isdb_df['pair_similarity'] = np.nan
    isdb_df['mps'] = np.nan
    isdb_df['pp'] = np.nan
    isdb_df['mps_pp'] = np.nan # sim最大，不一定是pp最大，所以要单独加一列
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
        isdb_df.loc[j, 'mps_pp'] = max_sim_pair[0]
        isdb_df.loc[j, 'pp'] = max_sim_pair[2]
    '''<Check> 输出+查看columns'''
    isdb_df.to_csv(isdb_file,index=None)

    '''计算1024_ecfp4'''
    # for i in trange(len(isdb_df)):
    #     try:
    #         edb_smile = isdb_df.loc[i, 'match_smiles']
    #         edb_mol = Chem.MolFromSmiles(edb_smile)
    #         smile = isdb_df.loc[i, 'smiles']
    #         mol = Chem.MolFromSmiles(smile)
    #         gnps_fp = AllChem.GetMorganFingerprintAsBitVect(edb_mol, radius=2, nBits=1024)
    #         fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=1024)
    #         sim = DataStructs.TanimotoSimilarity(gnps_fp, fp)
    #         isdb_df.loc[i, 'tanimoto'] = sim
    #     except:
    #         pass
    '''<Check>'''
    # isdb_df.to_csv('test.csv')

    '''Get identity in topK'''
    # top1, top2, top3, top4, top5, top10, top20 = 0, 0, 0, 0, 0, 0, 0
    # isdb_scans = list(Counter(isdb_df['row ID']))
    # index_isdb_to_keep1, index_isdb_to_keep5, index_isdb_to_keep20 = [], [], []
    # for isdb_scan in isdb_scans:
    #     isdb_df_slice = isdb_df[isdb_df['row ID'] == isdb_scan]
    #     isdb_df_slice_sorted = isdb_df_slice.sort_values('pp',ascending = False)
    # #     '''统计topK的identity数量'''
    #     top_1 = isdb_df_slice_sorted.head(1)['tanimoto'].values.tolist()
    #     top_5 = isdb_df_slice_sorted.head(5)['tanimoto'].values.tolist()
    #     top_20 = isdb_df_slice_sorted.head(20)['tanimoto'].values.tolist()
    #     if 1.0 in top_1:
    #         top1 += 1
    #     if 1.0 in top_5:
    #         top5 += 1
    #     if 1.0 in top_20:
    #         top20 += 1
    #
    #     index_top1 = isdb_df_slice_sorted.head(1).index
    #     index_isdb_to_keep1.extend(index_top1)
    #
    #     index_top5 = isdb_df_slice_sorted.head(5).index
    #     index_isdb_to_keep5.extend(index_top5)
    #
    #     index_top20 = isdb_df_slice_sorted.head(20).index
    #     index_isdb_to_keep20.extend(index_top20)
    #
    # print(top1, top5, top20)
    #
    # index_edb_to_keep1 = list(set(index_isdb_to_keep1))
    # index_edb_to_keep5 = list(set(index_isdb_to_keep5))
    # index_edb_to_keep20 = list(set(index_isdb_to_keep20))
    # edb_df_top1 = isdb_df.loc[index_edb_to_keep1]
    # edb_df_top5 = isdb_df.loc[index_edb_to_keep5]
    # edb_df_top20 = isdb_df.loc[index_edb_to_keep20]
    # edb_df_top1.to_csv('isdb_pp_top1.csv')
    # edb_df_top5.to_csv('isdb_pp_top5.csv')
    # edb_df_top20.to_csv('isdb_pp_top20.csv')





    print(f'Finished in {(time.time() - t) / 60:.2f} min')
