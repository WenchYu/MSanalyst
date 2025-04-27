
'''
The spectral data can be downloaded as Demo data 1 at (https://msanalyst.net/i/manuals)
The library searching can be repeated at (https://msanalyst.net/m/Molecular) using the above spectral data.
Filtering the GNPS library searching results based on the thresholds (matched peaks: 5; MS2 similarity: 0.1-0.9)

Input:
1. LibSearch_top100_0.1.csv
    The results can be downloaded at (url=https://gnps.ucsd.edu/ProteoSAFe/status.jsp?task=efe52ca63d664c06aab703e5d4f0bb64)
    Original name of the result file is MOLECULAR-LIBRARYSEARCH-V2-efe52ca6-view_all_annotations_DB-main copy.tsv
    Based on the original file, we manually add missing smiles and calculate the max common substructure similarity between standard compounds and matched smiles

2. std_info.csv
    Information about the standards

Output:
    ./data/GNPS/**.csv (Results filtered by similarity thresholds)

'''
import os,time

import numpy as np
import pandas as pd
from collections import Counter
from rdkit import Chem
from my_packages import cheminfo_tools
from tqdm import tqdm,trange

if __name__ == '__main__':
    t = time.time()

    output = './data/GNPS/'
    if not os.path.exists(output):
        os.mkdir(output)

    ''' Loading GNPS library searching result '''
    file = './data/LibSearch_top100_0.1.csv'
    df = pd.read_csv(file) # columns: ['#Scan#', 'Smiles']
    scans = list(Counter(df['#Scan#'])) # Get feature IDs

    std_df = pd.read_csv('./data/std_info.csv') # columns: ['name', 'smiles']

    ''' MCS calculating '''
    df['mcs'] = np.nan
    for i in trange(len(df)):
        feature_id = df.loc[i,'#Scan#']
        smile1 = df.loc[i,'Smiles'] # Smiles of reference spectra
        mol1 = Chem.MolFromSmiles(smile1)

        temp_std_df = std_df[std_df['name'] == feature_id]
        smile2 = temp_std_df['smiles'].values[0] # Smiles of standards
        mol2 = Chem.MolFromSmiles(smile2)

        mcs_similarity = cheminfo_tools.MCS(mol1,mol2)
        print(mcs_similarity)
        df.loc[i,'mcs1111'] = mcs_similarity

    df.to_csv(file,index=None) # save dataframe

    # Filtering by thresholds
    thresholds = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9]
    mps = [5]
    for mp in mps:
        for threshold in thresholds:
            df_filter = df[(df['MQScore'] >= threshold) & (df['SharedPeaks'] >= mp)]
            scans = list(Counter(df_filter['#Scan#']))
            index_to_keep = []
            for i in tqdm(scans,total= len(scans)):
                slice_df = df_filter[df_filter['#Scan#'] == i]
                max_index = slice_df['MQScore'].idxmax()
                index_to_keep.append(max_index)
            df_filtered = df.loc[index_to_keep].reset_index(drop=True)
            df_filtered.to_csv(f'{output}/{threshold}_{mp}.csv',index=None)

    print(f'Finished in {(time.time() - t) / 60:.2f} min')
