# -*- coding: utf-8 -*-
# @Time :2023/8/30 15:53
# @Auther :Yuwenchao
# @Software : PyCharm
'''
手动merge，计算出 MCS 的ISDB和EDB的计算结果成,merge.csv
1. 没有 mps 约束，sim, pp, mps 最大时的 MCS 分布, 均值或中位数
    最大值可能会有多个？直接取平均值或者取其中最大值

2. mps 约束时, sim, pp, sim + pp 最大时的 MCS 分布, 均值或中位数
    最大值可能会有多个？直接取平均值或者取其中最大值

3. 计算均值或中位数，绘制 density histogram

4. 画一个 3D 图，增加一个维度，表示 > 0.7 的数量

5.

'''
import os
import time

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from collections import Counter


if __name__ == '__main__':
    t = time.time()
    os.chdir('/Users/hehe/Desktop/HTS/MNA/#Fig2A-2F/match_results_mna_10ppm')

    file = '#merge.csv'
    df = pd.read_csv(file,index_col=0)


    '''
    1. 没有 mps 约束时
    '''
    # print(df.columns)
    # ids = list(Counter(df['row ID'])) # 获取 99 个feature ID 用于筛选
    #
    # # # demo用于写出测试循环中的合理代码
    # # temp = df[df['row ID']==ids[6]]
    # # max_sim = temp['pair_similarity'].max()
    # # max_sim_indices = temp[temp['pair_similarity'] == max_sim].index.tolist()
    # # print(np.mean(temp.loc[max_sim_indices, 'mcs']))
    #
    # mcs = [] # sim 或 pp 相同效果
    # mcs_mix = [] # sim + pp 同时match 的效果
    #
    # for id in ids:
    #     try:
    #         temp = df[df['row ID']==id]  # 无 mps约束的 temp_df
    #         max_sim = temp['pair_similarity'].max()
    #         max_sim_index = temp[temp['pair_similarity'] == max_sim].index.tolist()
    #
    #         max_pp = temp['pp'].max()
    #         max_pp_index = temp[temp['pp'] == max_pp].index.tolist()
    #
    #         max_mps = temp['mps'].max()
    #         max_mps_index = temp[temp['mps'] == max_mps].index.tolist()
    #
    #         sim_max = np.max(temp.loc[max_sim_index, 'mcs'])
    #         pp_max = np.max(temp.loc[max_pp_index, 'mcs'])
    #         mps_max = np.max(temp.loc[max_mps_index, 'mcs'])
    #
    #         # mcs.append(sim_max)
    #         # mcs.append(pp_max)
    #         # mcs.append(mps_max)
    #         mcs_mix.append(max([sim_max, pp_max]))
    #     except:
    #         pass

    '''
    2. mps 约束时
    '''
    print(df.columns)
    ids = list(Counter(df['row ID']))  # 获取 99 个feature ID 用于筛选

    # demo用于写出测试循环中的合理代码
    # temp = df[df['row ID']==ids[6]]
    # max_sim = temp['pair_similarity'].max()
    # max_sim_indices = temp[temp['pair_similarity'] == max_sim].index.tolist()
    # print(np.mean(temp.loc[max_sim_indices, 'mcs']))


    a,b,c,d = [],[],[],[]

    sim = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
    mps = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10 ]

    for i in mps:
        for j in sim:
            ann_count = 0
            mcs = []  # sim 或 pp 相同效果
            mcs_mix = []  # sim + pp 同时match 的效果
            for id in ids:
                try:
                    temp = df[(df['row ID'] == id) & (df['mps'] >= i) & ((df['pair_similarity'] >= j) | (df['pp'] >= j))]  # mps约束的 temp_df
                    # temp = temp[(temp['pair_similarity'] >= 0.7) | (temp['pp'] >= 0.7)]

                    max_sim = temp['pair_similarity'].max()
                    max_sim_index = temp[temp['pair_similarity'] == max_sim].index.tolist()

                    max_pp = temp['pp'].max()
                    max_pp_index = temp[temp['pp'] == max_pp].index.tolist()

                    max_mps = temp['mps'].max()
                    max_mps_index = temp[temp['mps'] == max_mps].index.tolist()

                    sim_max = np.max(temp.loc[max_sim_index, 'mcs'])
                    pp_max = np.max(temp.loc[max_pp_index, 'mcs'])
                    mps_max = np.max(temp.loc[max_mps_index, 'mcs'])

                    # mcs.append(sim_max)
                    mcs.append(pp_max)
                    # mcs.append(mps_max)
                    mcs_mix.append(max([sim_max, pp_max]))

                    if len(temp)>0:ann_count += 1
                except:
                    pass
            # print(f'sim={i} mps={j} ann_count={ann_count} mean={round(np.nanmean(mcs_mix),3)}')
            data = mcs_mix
            hit = sum(i >= 0.7  for i in data)
            per = round(hit/ann_count, 3)
            mean = round(np.nanmean(data), 3)
            print(f'{i}, {j}, {ann_count}, {hit}, {per}, {mean}')


    '''
    3. 
    '''
    # print(ann_count)
    # data = mcs_mix
    #
    # mean_value = round(np.nanmean(data),3)
    # median_value = round(np.nanmedian(data),3)
    # print(f'平均值：{mean_value}', f'中位数: {median_value}')
    # mcs2ser = pd.Series(data)
    # mu = np.mean(mcs2ser)
    # sigma = np.std(mcs2ser)
    # num_bins = 4
    # plt.figure(dpi=200)
    # plt.hist(mcs2ser
    #          , bins = num_bins
    #          , density = True
    #          , edgecolor = None
    #          , histtype = 'stepfilled'
    #          , color = 'steelblue'
    #          , alpha = 0.5
    #          , label = 'test'
    #          )
    # sns.kdeplot(mcs2ser)
    # ax = plt.gca()
    # # plt.legend()
    # plt.xlim(0, 1.5)  # 用你的范围替换 x_min 和 x_max
    # plt.ylim(0, 3)  # 用你的范围替换 y_min 和 y_max
    # plt.title(f'Mean: {mean_value}, Median: {median_value}')
    # plt.ylabel('Density')
    # plt.xlabel('MCS')
    # plt.show()
    # plt.savefig('',format = 'pdf')

    print(f'Finished in {(time.time() - t) / 60:.2f} min')
