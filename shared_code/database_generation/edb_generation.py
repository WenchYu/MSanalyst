
'''
Experimental library was collect from 'https://external.gnps2.org/gnpslibrary/ALL_GNPS.json'
This script was used for preprocess ALL_GNPS.json. And spectra acquired in positive mode was kept.
Demo data was provided (./data/edb_demo.json)
1. Use json to view the composition of the json file. {key : value}
2. Extract all "spectrum_id" and "peaks_json" to form the key-value pairs of {"spectrum_id" : "peaks_json"}
3. Generate 'edb_info.json' and 'edbMS1.csv', use "spectrum_id" as key.

'''
import os
import json
import time

import pandas as pd
from tqdm import tqdm,trange

if __name__ == '__main__':
    t = time.time()

    '''CFM-ID .log file preprocessing'''
    json_file = './data/edb_demo.json'
    with open(json_file, "r") as f:
        library_info = json.load(f)

    '''EDB MS2 library'''
    dict={}
    for i in trange(len(library_info)):
        dict[f'{library_info[i]["spectrum_id"]}'] ={
                                                    'pepmass' : f'{library_info[i]["Precursor_MZ"]}'
                                                    ,'charge' : f'{library_info[i]["Charge"]}'
                                                    ,'ion_mode' : f'{library_info[i]["Ion_Mode"]}'
                                                    ,'splash' : f'{library_info[i]["splash"]}'
                                                    ,'smiles': f'{library_info[i]["Smiles"]}'
                                                    ,'ms2': f'{library_info[i]["peaks_json"]}'
                                                    }

    with open ('./data/edb_info.json','w') as f:
        json.dump(dict,f)

    '''EDB MS1 library'''
    id,smiles,pepmass = [],[],[]
    for i in trange(len(library_info)):
        id.append(f'{library_info[i]["spectrum_id"]}')
        smiles.append(f'{library_info[i]["Smiles"]}')
        pepmass.append(f'{library_info[i]["Precursor_MZ"]}')
    data = {'id' : id, 'pepmass':pepmass, 'smiles' : smiles}
    df = pd.DataFrame(data)
    df.to_csv('./data/edbMS1.csv')


    print(f'Finish in {(time.time() - t) / 60:.2f}min')