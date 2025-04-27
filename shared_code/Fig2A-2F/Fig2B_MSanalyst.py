
'''
Fig 2B FBMN TopK annotations
input:
E_IS_merge.csv (Merged by IS_MS1match_std_quant.csv and E_MS1match_std_quant.csv)

output:
    Data for visualizing in prism
'''
import os
import time

import pandas as pd
from collections import Counter

if __name__ == '__main__':
    t = time.time()

    # Combining the 'E_MS1match_std_quant.csv' and 'IS_MS1match_std_quant.csv' and loading
    merge_file = './data/E_IS_merge.csv'
    merge_df = pd.read_csv(merge_file)
    print(merge_df.columns)

    # Counting the number of high similar pairs (MCS≥0.7) under topk spectral matches
    ids = list(Counter(merge_df['row ID'])) # Get a list of all features (45)
    top1, top2, top3, top4, top5, top10, top20 = 0, 0, 0, 0, 0, 0, 0
    index_to_keep1, index_to_keep5, index_to_keep20 = [], [], []
    for id in ids:
        df_slice = merge_df[merge_df['row ID'] == id]
        df_slice_sorted = df_slice.sort_values('pair_similarity', ascending=False) # Sort in descending order according to spectral similarity
        top_1 = df_slice_sorted.head(1)['mcs'].values.tolist()
        top_5 = df_slice_sorted.head(5)['mcs'].values.tolist()
        top_20 = df_slice_sorted.head(20)['mcs'].values.tolist()
        top1 += sum(i >= 0.7 for i in top_1)
        top5 += sum(i >= 0.7 for i in top_5)
        top20 += sum(i >= 0.7 for i in top_20)

        index_top1 = df_slice_sorted.head(1).index
        index_to_keep1.extend(index_top1)
        index_top5 = df_slice_sorted.head(5).index
        index_to_keep5.extend(index_top5)
        index_top20 = df_slice_sorted.head(20).index
        index_to_keep20.extend(index_top20)

    print(top1, top5, top20) # Data for visualization (figure 2B)

    # # Output the detailed infomation of topK results
    # index_edb_to_keep1 = list(set(index_to_keep1))
    # index_edb_to_keep5 = list(set(index_to_keep5))
    # index_edb_to_keep20 = list(set(index_to_keep20))
    # df_top1 = merge_df.loc[index_edb_to_keep1]
    # df_top5 = merge_df.loc[index_edb_to_keep5]
    # df_top20 = merge_df.loc[index_edb_to_keep20]
    # df_top1.to_csv('mix_top1.csv', index=None)
    # df_top5.to_csv('mix_top5.csv', index=None)
    # df_top20.to_csv('mix_top20.csv', index=None)
    print(f'Finished in {(time.time() - t) / 60:.2f} min')
