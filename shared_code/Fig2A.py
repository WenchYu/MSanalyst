# -*- coding: utf-8 -*-
# @Time :2023/6/18 16:42
# @Auther :Yuwenchao
# @Software : PyCharm
'''
Accumulation curve of GNPS and MNA results
'''
import time
import pandas as pd
from collections import Counter


if __name__ == '__main__':
    t = time.time()
    '''Data importing'''
    file = '/data/LibSearch_top100_0.1.tsv'
    df = pd.read_csv(file, sep='\t', encoding='gbk') # GNPS

    edb_file = '/data/edbMS1match_std_quant.csv' # MNA
    isdb_file = '/data/npMS1match_std_quant.csv'
    edb_df = pd.read_csv(edb_file)
    isdb_df = pd.read_csv(isdb_file)

    '''Screen and merge results of MNA from edb and isdb'''
    edb_df = edb_df[(edb_df['pair_similarity'] >= 0.1) | (edb_df['pp'] >= 0.1) & (edb_df['mps'] >= 1)]
    isdb_df = isdb_df[(isdb_df['pair_similarity'] >= 0.1) | (isdb_df['pp'] >= 0.1) & (isdb_df['mps'] >= 1)]
    columns_to_merge = ['row ID']
    merged_data = pd.concat(
        [isdb_df[columns_to_merge], edb_df[columns_to_merge]], axis=0)


    '''Counts'''
    scans = list(Counter(df['#Scan#']))
    counts = []
    for scan in scans:
        count = len(df[df['#Scan#']==scan])
        counts.append(count)
    mna_scans  = list(Counter(merged_data['row ID']))
    mna_counts = []
    for mna_scan in mna_scans:
        mna_count = len(merged_data[merged_data['row ID'] == mna_scan])
        mna_counts.append(mna_count)

    counts.sort()
    mna_counts.sort()
    y = []
    sum = 0
    for i in counts:
        sum+=i
        y.append(sum)

    y1 = []
    sum1 = 0
    for i in mna_counts:
        sum1 += i
        y1.append(sum1)

    print(y) # Y-axis coordinates for GNPS
    print(y1) # Y-axis coordinates for MNA

    print(f'Finished in {(time.time() - t) / 60:.2f} min')
