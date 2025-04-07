# -*- coding: utf-8 -*-
# @Time :2023/6/18 13:57
# @Auther :Yuwenchao
# @Software : PyCharm
'''
step-bar for GNPS
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


if __name__ == '__main__':
    t = time.time()
    '''data importing'''
    file = '/data/GNPSmerge.csv' # GNPS-FBMN matching results
    df = pd.read_csv(file,encoding='gbk')
    print(df.columns)

    '''TopK couning（MCS>0.7）'''
    ids = list(Counter(df['#Scan#']))
    top1, top2, top3, top4, top5, top10, top20 = 0, 0, 0, 0, 0, 0, 0
    index_to_keep1, index_to_keep5, index_to_keep20 = [], [], []
    for id in ids:
        df_slice = df[df['#Scan#'] == id]
        df_slice_sorted = df_slice.sort_values('MQScore', ascending=False) #descending
        top_1 = df_slice_sorted.head(1)['mcs'].values.tolist()
        top_5 = df_slice_sorted.head(5)['mcs'].values.tolist()
        top_20 = df_slice_sorted.head(20)['mcs'].values.tolist()
        top1 += sum(i >= 0.7 for i in top_1)
        top5 += sum(i >= 0.7 for i in top_5)
        top20 += sum(i >= 0.7 for i in top_20)

        index_top1 = df_slice_sorted.head(1).index #pick index of topK for specific investigation
        index_to_keep1.extend(index_top1)
        index_top5 = df_slice_sorted.head(5).index
        index_to_keep5.extend(index_top5)
        index_top20 = df_slice_sorted.head(20).index
        index_to_keep20.extend(index_top20)

    print(top1, top5, top20)  # Data for step-bar


    print(f'Finished in {(time.time() - t) / 60:.2f} min')
