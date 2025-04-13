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
    os.chdir('/Users/hehe/Desktop/HTS/MNA/#std_mix/match_results_mna_10ppm/mps5')

    threshold = 0.8
    mps = 5
    isdb_pp_file = f'npMS1match_std_quant.csv_max_pp_{threshold}_{mps}.csv'
    isdb_mq_file = f'npMS1match_std_quant.csv_max_sim_{threshold}_{mps}.csv'
    edb_pp_file = f'edbMS1match_std_quant.csv_max_pp_{threshold}_{mps}.csv'
    edb_mq_file = f'edbMS1match_std_quant.csv_max_sim_{threshold}_{mps}.csv'

    isdb_mq_df = pd.read_csv(isdb_mq_file)
    isdb_pp_df = pd.read_csv(isdb_pp_file)

    edb_mq_df = pd.read_csv(edb_mq_file)
    edb_pp_df = pd.read_csv(edb_pp_file)

    isdb_mq_df = isdb_mq_df.dropna(subset=['mcs'])
    isdb_pp_df = isdb_pp_df.dropna(subset=['mcs'])

    edb_mq_df = edb_mq_df.dropna(subset=['mcs'])
    edb_pp_df = edb_pp_df.dropna(subset=['mcs'])

    print('isdb sim', 'isdb pp', 'edb sim', 'edb pp')
    # print(f'{statistics.median(isdb_mq_df["mcs"])}\t{statistics.median(isdb_pp_df["mcs"])}\t{statistics.median(edb_mq_df["mcs"])}\t{statistics.median(edb_pp_df["mcs"])}')
    # print(
    #     f'{statistics.median(isdb_mq_df["mcs"])}')
    # print(
    #     f'{statistics.median(isdb_pp_df["mcs"])}')
    # print(
    #     f'{statistics.median(edb_mq_df["mcs"])}')
    # print(
    #     f'{statistics.median(edb_pp_df["mcs"])}')

    print(f'{np.average(isdb_mq_df["mcs"])}\t{np.average(isdb_pp_df["mcs"])}\t{np.average(edb_mq_df["mcs"])}\t{np.average(edb_pp_df["mcs"])}')





    print(f'Finished in {(time.time() - t) / 60:.2f} min')
