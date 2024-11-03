# -*- coding: utf-8 -*-
# @Time :2023/6/8 21:10
# @Auther :Yuwenchao
# @Software : PyCharm
'''
Because MNA additionally uses peak percentage and in-silico library,
the top 1 most similar result here is derived from the intersection of experimental and in-silico matches using peak percentage and modified cosine
Calculating the maximum common structure (MCS) similarity between structure of query features and matches
(High Similarity (>70%),Medium Similarity (45%-70%),Low Similarity (35%-45%), Not Similar (<35%), Not compared (under thresholds))

All match results were filtered based on the MCS,
retaining and counting only the top1 spectral similar match for each feature
'''
import os
import time

import numpy as np
import pandas as pd
from collections import Counter

if __name__ == '__main__':
    t = time.time()

    thresholds = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
    mps = 5

    os.chdir(f'/Users/hehe/Desktop/HTS/MNA/#std_mix/match_results_mna_10ppm/mps{mps}')

    for threshold in thresholds:
        isdb_pp_file = f'npMS1match_std_quant.csv_max_pp_{threshold}_{mps}.csv'
        isdb_mq_file = f'npMS1match_std_quant.csv_max_sim_{threshold}_{mps}.csv'
        edb_pp_file = f'edbMS1match_std_quant.csv_max_pp_{threshold}_{mps}.csv'
        edb_mq_file = f'edbMS1match_std_quant.csv_max_sim_{threshold}_{mps}.csv'

        isdb_pp_df = pd.read_csv(isdb_pp_file)
        isdb_mq_df = pd.read_csv(isdb_mq_file)
        edb_pp_df = pd.read_csv(edb_pp_file)
        edb_mq_df = pd.read_csv(edb_mq_file)

        columns_to_merge = ['row ID','mcs']
        merged_data = pd.concat(
            [isdb_pp_df[columns_to_merge], isdb_mq_df[columns_to_merge], edb_pp_df[columns_to_merge], edb_mq_df[columns_to_merge]], axis=0)
        # merged_data.to_csv('test.csv',index = None)

        scans = list(Counter(merged_data['row ID']))
        index_to_keep =[]
        for i in scans:
            slice_df = merged_data[merged_data['row ID'] == i]
            index_max = slice_df['mcs'].idxmax()
            if not np.isnan(index_max):
                index_to_keep.append(index_max)
        df_filtered = merged_data.loc[index_to_keep].drop_duplicates(subset='row ID').reset_index(drop=True)
        # df_filtered.to_csv(f'test111.csv', index=None)

        top1, H, M, L, NS, NC = 0, 0, 0, 0, 0, 0
        for i in range(len(df_filtered)):
            mcs = df_filtered.loc[i, 'mcs']
            if mcs == 1:
                top1 += 1
            elif mcs >= 0.7 and mcs < 1:
                H += 1
            elif mcs >= 0.45 and mcs < 0.7:
                M += 1
            elif mcs >= 0.35 and mcs < 0.45:
                L += 1
            else:
                NS += 1
        NC = 99 - H - M - L - NS - top1
        print(f'{NC}\t{NS}\t{L}\t{M}\t{H}\t{top1}')


    print(f'Finished in {(time.time() - t) / 60:.2f} min')
