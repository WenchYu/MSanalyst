# -*- coding: utf-8 -*-
# @Time :2023/6/7 23:03
# @Auther :Yuwenchao
# @Software : PyCharm
'''
Based on the MNA results,
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
from tqdm import trange

def MSC(mol1, mol2):
    # Define parameters using variables
    # min_num_atoms = 4
    # threshold = 0.5
    # atom_compare = rdFMCS.AtomCompare.CompareAny
    # bond_compare = rdFMCS.BondCompare.CompareAny
    #
    # # Use the parameters in your code
    # mcs_params = rdFMCS.MCSParameters()
    # mcs_params.minNumAtoms = min_num_atoms
    # mcs_params.threshold = threshold
    # mcs_params.atomCompare = atom_compare
    # mcs_params.bondCompare = bond_compare
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


    '''pp 或 MQ file'''
    mps = 3
    thresholds = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
    os.chdir(f'/Users/hehe/Desktop/HTS/MNA/#std_mix/match_results_mna_10ppm/mps{mps}')

    for threshold in thresholds:
        mna_file = f'edbMS1match_std_quant.csv_max_pp_{threshold}_{mps}.csv'
        mna_df = pd.read_csv(mna_file)

        # mna_df['mcs'] = np.nan
        # for i in trange(len(mna_df)):
        #     try:
        #         smile1 = mna_df.loc[i, 'smiles']
        #         mol1 = Chem.MolFromSmiles(smile1)
        #         smile2 = mna_df.loc[i, 'match_smiles']
        #         mol2 = Chem.MolFromSmiles(smile2)
        #         mna_df.loc[i, 'mcs'] = MSC(mol1, mol2)
        #     except:
        #         pass
        # mna_df.to_csv(mna_file,index=None)

        top1, H, M, L, NS, NC = 0, 0, 0, 0, 0, 0
        for i in range(len(mna_df)):
            mcs = mna_df.loc[i, 'mcs']
            if mcs == 1:
                top1 += 1
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
