
'''
Fig 2B GNPS TopK annotations
input:
LibSearch_top100_0.1.csv (GNPS library searching results)

output:
    Data for visualizing in prism



'''
import os
import time
import pandas as pd
import numpy as np
from tqdm import trange
from collections import Counter


if __name__ == '__main__':
    t = time.time()

    # Loading GNPS library searching result
    file = './data/LibSearch_top100_0.1.csv'
    df = pd.read_csv(file)

    # Counting the number of high similar pairs under topk (MCS≥0.7)
    ids = list(Counter(df['#Scan#']))  # Get a list of all features (45)
    top1, top2, top3, top4, top5, top10, top20 = 0, 0, 0, 0, 0, 0, 0
    index_to_keep1, index_to_keep5, index_to_keep20 = [], [], []
    for id in ids:
        df_slice = df[df['#Scan#'] == id]
        df_slice_sorted = df_slice.sort_values('MQScore', ascending=False) # Sort in descending order according to spectral similarity

        top_1 = df_slice_sorted.head(1)['mcs'].values.tolist()
        top_5 = df_slice_sorted.head(5)['mcs'].values.tolist()
        top_20 = df_slice_sorted.head(20)['mcs'].values.tolist()
        top1 += sum(i >= 0.7 for i in top_1)
        top5 += sum(i >= 0.7 for i in top_5)
        top20 += sum(i >= 0.7 for i in top_20)

        # pick index of topK for manual investigation
        index_top1 = df_slice_sorted.head(1).index
        index_to_keep1.extend(index_top1)
        index_top5 = df_slice_sorted.head(5).index
        index_to_keep5.extend(index_top5)
        index_top20 = df_slice_sorted.head(20).index
        index_to_keep20.extend(index_top20)

    print(top1, top5, top20)  # Data for visualization (figure 2B)

    # # Output the detailed infomation of topK results
    # index_edb_to_keep1 = list(set(index_to_keep1))
    # index_edb_to_keep5 = list(set(index_to_keep5))
    # index_edb_to_keep20 = list(set(index_to_keep20))
    # df_top1 = df.loc[index_edb_to_keep1]
    # df_top5 = df.loc[index_edb_to_keep5]
    # df_top20 = df.loc[index_edb_to_keep20]
    # df_top1.to_csv('fbmn_top1.csv', index=None)
    # df_top5.to_csv('fbmn_top5.csv', index=None)
    # df_top20.to_csv('fbmn_top20.csv', index=None)

    print(f'Finished in {(time.time() - t) / 60:.2f} min')
