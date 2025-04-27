
'''
The spectral data can be downloaded as Demo data 1 at (https://msanalyst.net/i/manuals)
The library searching can be repeated at (https://msanalyst.net/m/Molecular) using the above spectral data.
Filtering the MSanalyst results based on the thresholds (matched peaks: 5; MS2 similarity/peak percentage: 0.1-0.9)

input:
    E_MS1match_std_quant.csv
    IS_MS1match_std_quant.csv

output:
    ./data/MSanalyst/E_**_pp.csv (Experimental library search results + peak percentage filtering)
    ./data/MSanalyst/E_**_sim.csv (Experimental library search results + modified_cosine filtering)
    ./data/MSanalyst/IS_**_pp.csv (In-silico library search results + peak percentage filtering)
    ./data/MSanalyst/IS_**_sim.csv (In-silico library search results + modified_cosine filtering)

'''

import os
import time

import numpy as np
import pandas as pd

from tqdm import tqdm,trange
from collections import Counter
if __name__ == '__main__':
    t = time.time()

    output = './data/MSanalyst/'
    if not os.path.exists(output):
        os.mkdir(output)

    '''Loading msanalyst searching result (experimental library)'''
    mq_file = './data/E_MS1match_std_quant.csv'
    mq_df = pd.read_csv(mq_file) #
    type = 'max_pp' # switch to 'max_sim' and run it again
    thresholds = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9]
    mps = 5 # 3，4，5，6

    # Filtering the results by peak percentage and modified_cosine respectively
    for threshold in thresholds:
        mq_df_filtered = mq_df[(mq_df['pp'] >= threshold) & (mq_df['mps'] >= mps)]
        # mq_df_filtered = mq_df[(mq_df['pair_similarity']>=threshold) & (mq_df['mps']>=mps)] #
        scans = list(Counter(mq_df_filtered['row ID']))
        index_to_keep = []
        for i in tqdm(scans, total = len(scans)):
            feature_id = i
            slice_df = mq_df_filtered[mq_df['row ID'] == feature_id]
            max_index = slice_df['pp'].idxmax()
            # max_index = slice_df['pair_similarity'].idxmax()
            if not np.isnan(max_index):
                index_to_keep.append(max_index)
        df_filtered = mq_df.loc[index_to_keep].reset_index(drop=True)
        df_filtered.to_csv(f'{output}{os.path.basename(mq_file)}_{type}_{threshold}_{mps}.csv',index=None)

    '''Loading msanalyst searching result (in silico library)'''
    is_file = './data/IS_MS1match_std_quant.csv'
    is_df = pd.read_csv(is_file)
    is_type = 'max_sim' # types = 'max_pp'
    is_df[type],is_df['max_mps'] = np.nan, np.nan

    for threshold in thresholds:
        for j in trange(len(is_df)):
            pairs = [
                (is_df.loc[j, 'mps0'], is_df.loc[j, 'pair_similarity0'], is_df.loc[j, 'pp0']),
                (is_df.loc[j, 'mps1'], is_df.loc[j, 'pair_similarity1'], is_df.loc[j, 'pp1']),
                (is_df.loc[j, 'mps2'], is_df.loc[j, 'pair_similarity2'], is_df.loc[j, 'pp2'])]

            # Retrieve the results with maximum modified_cosine or peak_percentage
            # switch [1] for similarity and [2] for peak peak percentage

            # max_pair = max(pairs, key=lambda x: x[2]) # [0] mps [1] pair_similarity [2] pp
            max_pair = max(pairs, key=lambda x: x[1])  # [0] mps [1] pair_similarity [2] pp

            is_df.loc[j, 'max_mps'] = max_pair[0]
            # is_df.loc[j, type] = max_pair[2]
            is_df.loc[j, type] = max_pair[1]

        index_to_keep=[]
        np_df_filtered = is_df[(is_df[type] >= threshold) & (is_df['max_mps']>=mps) ]
        scans = list(Counter(np_df_filtered['row ID']))
        for i in tqdm(scans, total = len(scans)):
            slice_df = np_df_filtered[is_df['row ID']==i]
            max_index = slice_df[type].idxmax()
            if not np.isnan(max_index):
                index_to_keep.append(max_index)
        df_filtered = is_df.loc[index_to_keep].reset_index(drop=True)
        df_filtered.to_csv(f'{output}{os.path.basename(is_file)}_{is_type}_{threshold}_{mps}.csv',index=None)


    print(f'Finished in {(time.time() - t) / 60:.2f} min')
