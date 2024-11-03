# -*- coding: utf-8 -*-
# @Time :2023/6/7 20:14
# @Auther :Yuwenchao
# @Software : PyCharm
'''
Based on the GNPS library searching results,
calculating the maximum common structure (MCS) similarity between structure of query features and matches
(High Similarity (>70%),Medium Similarity (45%-70%),Low Similarity (35%-45%), Not Similar (<35%), Not compared (under thresholds))

All match results were filtered based on the MCS,
retaining and counting only the top1 spectral similar match for each feature
'''
import os
import time
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import rdFMCS

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

if __name__ == '__main__':
    t = time.time()
    mps = 5
    thresholds = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
    for threshold in thresholds:
        gnps_file = f'LibSearch_top100_0.1.tsv_{threshold}_{mps}.csv'
        gnps_df = pd.read_csv(gnps_file)
        gnps_df['mcs'] = np.nan
        for i in range(len(gnps_df)):
            try:
                smile1 = gnps_df.loc[i,'smiles']
                mol1 = Chem.MolFromSmiles(smile1)
                smile2 = gnps_df.loc[i,'Smiles']
                mol2 = Chem.MolFromSmiles(smile2)
                gnps_df.loc[i,'mcs'] =MSC(mol1,mol2)
            except:pass
        gnps_df.to_csv(gnps_file, index=None)

        top1, H, M, L, NS, NC = 0, 0, 0, 0, 0, 0
        for i in range(len(gnps_df)):
            mcs = gnps_df.loc[i, 'mcs']
            if mcs == 1:
                top1 +=1
            elif mcs >= 0.7 and mcs < 1:
                H += 1
            elif mcs >= 0.45 and mcs < 0.7:
                M += 1
            elif mcs >= 0.35 and mcs < 0.45:
                L += 1
            else:
                NS += 1
        NC = 99 - H - M - L - NS - top1
        print(f'{NC}\t{NS}\t{L}\t{M}\t{H}\t{top1}')

    print(f'Finished in {(time.time() - t) / 60:.2f} min')
