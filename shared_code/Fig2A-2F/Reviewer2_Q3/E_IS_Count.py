# -*- coding: utf-8 -*-
# @Time :2025/1/24 23:00
# @Auther :Yuwenchao
# @Software : PyCharm
'''
根据CCMSL，no-CCMSL 数量统计

'''
import time
import pandas as pd


def IdConts(df, mcs_threshold=0.7):
    """
    统计给定 DataFrame 中 match_id 列包含和不包含 'CCMS' 的行数量，基于指定的 mcs 阈值进行过滤。

    参数:
    df (pd.DataFrame): 输入的 DataFrame。
    mcs_threshold (float): 过滤 mcs 列的阈值，默认值为 0.7。

    返回:
    tuple: 包含 'CCMS' 的行数量 (edb_counts) 和不包含 'CCMS' 的行数量 (isdb_counts)。
    """
    # 过滤 mcs 列值大于等于阈值的行
    filtered_df = df[df['mcs'] > mcs_threshold]
    # filtered_df = filtered_df.drop_duplicates(subset=drop_duplicates_column)

    # 统计包含和不包含 'CCMS' 的行数
    edb_counts = filtered_df['match_id'].str.contains('CCMS').sum()
    isdb_counts = (~filtered_df['match_id'].str.contains('CCMS')).sum()

    return edb_counts, isdb_counts


if __name__ == '__main__':
    t = time.time()

    top1df = pd.read_csv('mix_top1.csv')
    top5df = pd.read_csv('mix_top5.csv')
    top20df = pd.read_csv('mix_top20.csv')

    top1eC, top1isC = IdConts(top1df,0.7)
    top5eC, top5isC = IdConts(top5df, 0.7)
    top20eC, top20isC = IdConts(top20df, 0.7)

    print(top1eC,top1isC)
    print(top5eC, top5isC)
    print(top20eC, top20isC)

    print(f'Finished in {(time.time() - t) / 60:.2f} min')
