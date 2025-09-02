
# -*- coding: utf-8 -*-
# @Time :2025/6/12 00:10
# @Auther :Yuwenchao
# @Software : PyCharm
'''

'''
import sys,os,time,json
sys.path.append('../')
from my_packages import functions_new, cheminfo_tools, evaluation,peaktools
import numpy as np
from joblib import Parallel, delayed
from tqdm import tqdm,trange

def compute_mcs(idx1, idx2, CCMSID1, CCMSID2, GNPS_INFO):
    SMILE1 = GNPS_INFO[CCMSID1]['CANONSMILES']
    SMILE2 = GNPS_INFO[CCMSID2]['CANONSMILES']
    MCS = cheminfo_tools.mcs(SMILE1, SMILE2)
    return idx1, idx2, MCS

if __name__ == '__main__':
    t = time.time()
    # Shared peak matrices
    GNPS_INFO = functions_new.json_load('../msdb/GNPSLIBRARY_250514/GNPS-LIBRARY-INFO.json')
    CCMSIDs = np.load('../msdb/data/idlist/H_qtof_non-redundant_CCMSIDs.npy')[:15]

    MCS_MATRIX = np.zeros((len(CCMSIDs), (len(CCMSIDs))))
    # 创建一个任务列表，包含所有需要计算的 (idx1, idx2) 对
    tasks = []
    for idx1 in range(len(CCMSIDs)):
        CCMSID1 = CCMSIDs[idx1]
        for idx2 in range(idx1):
            CCMSID2 = CCMSIDs[idx2]
            tasks.append(delayed(compute_mcs)(idx1, idx2, CCMSID1, CCMSID2, GNPS_INFO))

    results = Parallel(n_jobs=6, verbose=10)(tasks)

    for idx1, idx2, MCS in results:
        MCS_MATRIX[idx1, idx2] = MCS
        MCS_MATRIX[idx2, idx1] = MCS

    np.save('mcs_similarity_matrix.npy',MCS_MATRIX)

    # Save indices <str> : CCMSIDs
    INDEXtoCCMSID = {index: CCMSID for index, CCMSID in enumerate(CCMSIDs)}
    with open (f'mcs_INDEXtoCCMS.json', 'w') as f2:
        json.dump(INDEXtoCCMSID,f2)

    print(f'Finished in {(time.time() - t) / 60:.2f} min')
