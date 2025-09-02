# -*- coding: utf-8 -*-
# @Time :2025/8/1 17:34
# @Auther :Yuwenchao
# @Software : PyCharm
'''

'''
import time,os
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib
matplotlib.use('Agg')


if __name__ == '__main__':
    t = time.time()

    # 1D hist
    DIR = f'../msdb/data/hqtof/SpecSimMatrix/'
    ALGORITHMs = ['cosine']
    for ALGO in ALGORITHMs:

        SIMILARITY = []  # Hist data
        SPEC_SIM_MATRIX = np.load(os.path.join(DIR, f'{ALGO}.npy'))
        for idx1 in range(len(SPEC_SIM_MATRIX)):
            for idx2 in range(idx1):
                SIMILARITY.append(SPEC_SIM_MATRIX[idx1, idx2])
                
        plt.figure(figsize=(6, 4))  # 设置图形大小
        sns.kdeplot(SIMILARITY, fill=True, color='blue', label=ALGO)  # 绘制 KDE 图，填充颜色
        plt.title(f'KDE Plot of Similarity for {ALGO}')  # 设置标题
        plt.xlabel('Similarity')  # 设置 x 轴标签
        plt.ylabel('Density')  # 设置 y 轴标签
        plt.legend()  # 显示图例
        plt.show()  # 显示图形
    print(f'Finished in {(time.time() - t) / 60:.2f} min')
