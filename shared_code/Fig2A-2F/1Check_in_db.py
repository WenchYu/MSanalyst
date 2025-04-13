
'''
检查lxn给的化合物是否在gnps，coconut，npatls，cmnpd的数据库里
结构数据库部分手动进行
以下代码是查找gnps
'''

import os
import time

import numpy as np
import pandas as pd
from tqdm import tqdm,trange
from rdkit import Chem

if __name__ == '__main__':
    t = time.time()
    os.chdir('/Users/hehe/Desktop/HTS/msdb')
    db_file = 'gnps.csv'
    db_df = pd.read_csv(db_file)
    print(db_df.columns)


    '''将所有 smiles 转换成 canonical smile'''
    # db_df['canon'] = np.nan
    # for i in trange(len(db_df)):
    #     try:
    #         db_smile = db_df.loc[i,'smiles']
    #         db_mol = Chem.MolFromSmiles(db_smile)
    #         canon_smiles1 = Chem.MolToSmiles(db_mol, isomericSmiles=True, canonical=True)
    #         db_df.loc[i,'canon'] = canon_smiles1
    #     except:
    #         pass
    # db_df.to_csv('gnps.csv',index=None)

    file = 'test.csv'
    df = pd.read_csv(file)
    print(df.columns)

    # df['canon'] = np.nan
    # for i in trange(len(df)):
    #     try:
    #         smile = df.loc[i,'smiles']
    #         mol = Chem.MolFromSmiles(smile)
    #         canon_smiles = Chem.MolToSmiles(mol, isomericSmiles=True, canonical=True)
    #         df.loc[i,'canon'] = canon_smiles
    #     except:
    #         pass
    # df.to_csv('test1111.csv',index=None)

    df['gnps_id'] = np.nan
    for i in trange(len(df)):
        smile = df.loc[i,'canon']
        if smile in db_df['canon'].values:
            indices = db_df.id[db_df['canon'] == smile].values.tolist()
            print(indices[0])
            df.loc[i, 'gnps_id'] = indices[0]

    df.to_csv('test111.csv',index=None)

    print(f'Finished in {(time.time() - t) / 60:.2f} min')
