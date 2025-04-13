# -*- coding: utf-8 -*-
# @Time :2023/9/4 15:39
# @Auther :Yuwenchao
# @Software : PyCharm
'''

'''
import os
import time
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.cm import get_cmap
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

if __name__ == '__main__':
    t = time.time()
    os.chdir('/Users/hehe/Desktop/HTS/MNA/#Fig2A-2F/match_results_mna_10ppm')
    file = 'sim_pp.csv'
    df = pd.read_csv(file)
    print(df.columns)


    temp = df[df['mps'] == 1]
    group1_x = temp['sim']
    group1_y = temp['hits']
    group1_z = temp['mcs']

    temp = df[df['mps'] == 2]
    group2_x = temp['sim']
    group2_y = temp['hits']
    group2_z = temp['mcs']

    temp = df[df['mps'] == 3]
    group3_x = temp['sim']
    group3_y = temp['hits']
    group3_z = temp['mcs']

    temp = df[df['mps'] == 4]
    group4_x = temp['sim']
    group4_y = temp['hits']
    group4_z = temp['mcs']

    temp = df[df['mps'] == 5]
    group5_x = temp['sim']
    group5_y = temp['hits']
    group5_z = temp['mcs']

    temp = df[df['mps'] == 6]
    group6_x = temp['sim']
    group6_y = temp['hits']
    group6_z = temp['mcs']

    temp = df[df['mps'] == 7]
    group7_x = temp['sim']
    group7_y = temp['hits']
    group7_z = temp['mcs']

    temp = df[df['mps'] == 8]
    group8_x = temp['sim']
    group8_y = temp['hits']
    group8_z = temp['mcs']

    temp = df[df['mps'] == 9]
    group9_x = temp['sim']
    group9_y = temp['hits']
    group9_z = temp['mcs']

    temp = df[df['mps'] == 10]
    group10_x = temp['sim']
    group10_y = temp['hits']
    group10_z = temp['mcs']

    cmap = get_cmap('tab10', 10)

    # 创建一个图形和3D坐标轴
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    # 绘制第一组散点
    # ax.plot(group1_x, group1_y, group1_z, c=cmap(0), marker='o', label='mps=1', markersize=3, linewidth=0.5)
    # ax.plot(group2_x, group2_y, group2_z, c=cmap(1), marker='o', label='mps=2', markersize=3, linewidth=0.5)
    # ax.plot(group3_x, group3_y, group3_z, c=cmap(2), marker='o', label='mps=3', markersize=3, linewidth=0.5)
    # ax.plot(group4_x, group4_y, group4_z, c=cmap(3), marker='o', label='mps=4', markersize=3, linewidth=0.5)
    ax.plot(group5_x, group5_y, group5_z, c=cmap(4), marker='o', label='mps=5', markersize=3, linewidth=0.5)
    ax.plot(group6_x, group6_y, group6_z, c=cmap(5), marker='o', label='mps=6', markersize=3, linewidth=0.5)
    ax.plot(group7_x, group7_y, group7_z, c=cmap(6), marker='o', label='mps=7', markersize=3, linewidth=0.5)
    # ax.plot(group8_x, group8_y, group8_z, c=cmap(7), marker='o', label='mps=8', markersize=3, linewidth=0.5)
    # ax.plot(group9_x, group9_y, group9_z, c=cmap(8), marker='o', label='mps=9', markersize=3, linewidth=0.5)
    # ax.plot(group10_x, group10_y, group10_z, c=cmap(9), marker='o', label='mps=10', markersize=3, linewidth=0.5)




    # 设置坐标轴标签
    ax.set_xlabel('MS/MS similarity')
    ax.set_ylabel('Hits')
    # ax.set_ylim(10, 20)  # 这里将y轴限制在0到1之间
    ax.set_zlabel('MCS')
    ax.set_zlim(0.5, 0.9)

    # ax.grid(False) # 关闭网格线

    # 将图例移到更右更上
    ax.legend(loc='right', bbox_to_anchor=(1, 1),fontsize=5)

    ax.view_init(elev=20, azim=145) # elev仰角, azim 方位角

    # 添加图例
    # ax.legend()

    # 显示图形
    plt.show()

    print(f'Finished in {(time.time() - t) / 60:.2f} min')
