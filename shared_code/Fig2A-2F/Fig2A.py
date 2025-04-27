
'''
Figure 2A. Accumulation curve of GNPS and MSanalyst results

input:
LibSearch_top100_0.1.csv (GNPS library searching results)

E_MS1match_std_quant.csv (Experimental library searching results)
IS_MS1match_std_quant.csv (In-silico library searching results)

output:
    Data for visualizing in prism
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
    # Loading GNPS library searching result
    file = './data/LibSearch_top100_0.1.csv'
    df = pd.read_csv(file)

    # Loading MSanalyst library searching result
    edb_file = './data/E_MS1match_std_quant.csv'
    isdb_file = './data/IS_MS1match_std_quant.csv'
    edb_df = pd.read_csv(edb_file)
    isdb_df = pd.read_csv(isdb_file)

    # Filtering and merging the MS1 match result
    edb_df = edb_df[(edb_df['pair_similarity'] >= 0.1) | (edb_df['pp'] >= 0.1) & (edb_df['mps'] >= 1)]
    isdb_df = isdb_df[(isdb_df['pair_similarity'] >= 0.1) | (isdb_df['pp'] >= 0.1) & (isdb_df['mps'] >= 1)]
    columns_to_merge = ['row ID']
    merged_data = pd.concat(
        [isdb_df[columns_to_merge], edb_df[columns_to_merge]], axis=0)

    # GNPS MS1 match counts
    scans = list(Counter(df['#Scan#']))
    counts = []
    for scan in scans:
        count = len(df[df['#Scan#']==scan])
        counts.append(count)

    # MSanalyst MS1 match counts
    msanalyst_scans  = list(Counter(merged_data['row ID']))
    msanalyst_counts = []
    for mna_scan in msanalyst_scans:
        mna_count = len(merged_data[merged_data['row ID'] == mna_scan])
        msanalyst_counts.append(mna_count)

    counts.sort()
    msanalyst_counts.sort()
    y = []
    sum = 0
    for i in counts:
        sum+=i
        y.append(sum)

    y1 = []
    sum1 = 0
    for i in msanalyst_counts:
        sum1 += i
        y1.append(sum1)

    print(f'GNPS ms1 matched {len(y)} features')
    print(f'MSanalyst ms1 match{len(y1)} features')
    print(f'Data for GNPS accumulative curve {y}')
    print(f'Data for MSanalyst accumulative curve {y1}') #

    print(f'Finished in {(time.time() - t) / 60:.2f} min')
