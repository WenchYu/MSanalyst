# -*- coding: utf-8 -*-
# @Time :2023/6/11 21:58
# @Auther :Yuwenchao
# @Software : PyCharm
'''

'''
import os
import statistics
import time

import numpy as np
import pandas as pd

if __name__ == '__main__':
    t = time.time()
    os.chdir('/Users/hehe/Desktop/HTS/MNA/#std_mix/GNPS_LibSearch/mps6')

    threshold = 0.8
    mps = 6
    gnps_file = f'LibSearch_top100_0.1.tsv_{threshold}_{mps}.csv'

    gnps_df = pd.read_csv(gnps_file)
    gnps_df = gnps_df.dropna(subset=['mcs'])

    # print(f'{statistics.median(gnps_df["mcs"])}')
    print(f'{np.average(gnps_df["mcs"])}')


    print(f'Finished in {(time.time() - t) / 60:.2f} min')
