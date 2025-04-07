# -*- coding: utf-8 -*-
# @Time :2023/6/18 16:42
# @Auther :Yuwenchao
# @Software : PyCharm
'''
step-bar for GNPS
'''

import time
import pandas as pd
from collections import Counter

if __name__ == '__main__':
    t = time.time()
    '''Data importing'''
    merge_file = '/data/msanalystmerge.csv' # edb  and isdb match results
    merge_df = pd.read_csv(merge_file)
    print(merge_df.columns)



    '''TopK Counts（MCS>0.7）'''
    ids = list(Counter(merge_df['row ID']))
    top1, top2, top3, top4, top5, top10, top20 = 0, 0, 0, 0, 0, 0, 0
    index_to_keep1, index_to_keep5, index_to_keep20 = [], [], []
    for id in ids:
        df_slice = merge_df[merge_df['row ID'] == id]
        df_slice_sorted = df_slice.sort_values('pair_similarity', ascending=False) #descending
        top_1 = df_slice_sorted.head(1)['mcs'].values.tolist()
        top_5 = df_slice_sorted.head(5)['mcs'].values.tolist()
        top_20 = df_slice_sorted.head(20)['mcs'].values.tolist()
        top1 += sum(i > 0.7 for i in top_1)
        top5 += sum(i > 0.7 for i in top_5)
        top20 += sum(i > 0.7 for i in top_20)

        index_top1 = df_slice_sorted.head(1).index
        index_to_keep1.extend(index_top1)
        index_top5 = df_slice_sorted.head(5).index
        index_to_keep5.extend(index_top5)
        index_top20 = df_slice_sorted.head(20).index
        index_to_keep20.extend(index_top20)

    print(top1, top5, top20) # Data for step-bar

    print(f'Finished in {(time.time() - t) / 60:.2f} min')
