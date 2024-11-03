# -*- coding: utf-8 -*-
# @Time :2023/6/8 10:48
# @Auther :Yuwenchao
# @Software : PyCharm
'''
Dividing the GNPS library searching results(url=https://gnps.ucsd.edu/ProteoSAFe/status.jsp?task=efe52ca63d664c06aab703e5d4f0bb64)
based on the thresholds (MPC ranging from 3 to 6; MS2 similarity ranging from 0.1-0.9)
Output folder/files used for fig2C
'''

import time

import pandas as pd
from collections import Counter
from tqdm import tqdm
if __name__ == '__main__':
    t = time.time()

    file = '/data/LibSearch_top100_0.1.tsv'
    df = pd.read_csv(file,sep = '\t',encoding='gbk')
    scans = list(Counter(df['#Scan#']))
    print(len(scans))

    print(df.columns)
    thresholds = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9]
    mps = [5]

    for mp in mps:
        for threshold in thresholds:
            df_filter = df[(df['MQScore'] >= threshold) & (df['SharedPeaks'] >= mp)]
            scans = list(Counter(df_filter['#Scan#']))
            index_to_keep = []
            for i in tqdm(scans,total= len(scans)):
                slice_df = df_filter[df_filter['#Scan#'] == i]
                max_index = slice_df['MQScore'].idxmax()
                index_to_keep.append(max_index)
            df_filtered = df.loc[index_to_keep].reset_index(drop=True)
            df_filtered.to_csv(f'{file}_{threshold}_{mp}.csv',index=None)

    print(f'Finished in {(time.time() - t) / 60:.2f} min')
