# -*- coding: utf-8 -*-
# @Time :2025/1/24 23:40
# @Auther :Yuwenchao
# @Software : PyCharm
'''

'''
import time

import pandas as pd

if __name__ == '__main__':
    t = time.time()
    df = pd.read_csv('std_mix_node_info.csv',encoding='gbk')
    ACounts = df['level'].str.contains('A').sum()
    B1Counts = df['level'].str.contains('B1').sum()
    B2Counts = df['level'].str.contains('B2').sum()
    C1Counts = df['level'].str.contains('C1').sum()
    C2Counts = df['level'].str.contains('C2').sum()

    print(ACounts, B1Counts, B2Counts, C1Counts, C2Counts)

    print(f'Finished in {(time.time() - t) / 60:.2f} min')
