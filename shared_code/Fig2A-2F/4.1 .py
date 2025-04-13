# -*- coding: utf-8 -*-
# @Time :2023/8/30 14:01
# @Auther :Yuwenchao
# @Software : PyCharm
'''
calculate_similarity 会受空 molecule的影响，空smiles返回 0 即可

'''

import os
import time
import pandas as pd
import numpy as np
from rdkit import Chem
from tqdm import trange, tqdm
from rdkit.Chem import rdFMCS
from joblib import Parallel, delayed

def MSC(mol1, mol2):
    mcs = rdFMCS.FindMCS([mol1, mol2]
                         , bondCompare=rdFMCS.BondCompare.CompareOrder
                         , atomCompare=rdFMCS.AtomCompare.CompareAny
                         , maximizeBonds = False
                         , ringMatchesRingOnly=False
                         , matchValences=False
                         , timeout=10
                         )

    mcs_num_bonds = mcs.numBonds
    mol1_num_bonds = mol1.GetNumBonds()
    mol2_num_bonds = mol2.GetNumBonds()

    similarity = mcs_num_bonds / ((mol1_num_bonds + mol2_num_bonds) - mcs_num_bonds)
    return similarity

def calculate_similarity(index, row):
    smile = row['smiles']
    mol1 = Chem.MolFromSmiles(smile)

    match_smile = row['match_smiles']
    try:
        mol2 = Chem.MolFromSmiles(match_smile)
        similarity = MSC(mol1, mol2)
    except:
        similarity = 0

    return index, similarity


if __name__ == '__main__':
    t = time.time()
    os.chdir('/Users/hehe/Desktop/HTS/MNA/#Fig2A-2F/match_results_mna_10ppm')

    file = 'edbMS1match_test.csv'
    df = pd.read_csv(file)
    print(df.columns)

    # 获取CPU核心数量
    num_cores = os.cpu_count()

    # 使用joblib的Parallel并行处理任务
    results = Parallel(n_jobs=num_cores)(
        delayed(calculate_similarity)(i, row) for i, row in tqdm(df.iterrows(),total = len(df))
    )

    # 更新DataFrame中的计算结果
    for index, similarity in results:
        df.loc[index, 'mcs'] = similarity

    df.to_csv('test.csv')


print(f'Finished in {(time.time() - t) / 60:.2f} min')
