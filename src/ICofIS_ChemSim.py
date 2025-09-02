# -*- coding: utf-8 -*-
# @Time :2025/7/28 00:31
# @Auther :Yuwenchao
# @Software : PyCharm
'''
Calculate chemical similarity pairs in ISDB from each other to determine their Inner consistency
[FPsim2 documentation](https://chembl.github.io/FPSim2/user_guide/install/)
'''
import json
import os,time,argparse,sys,random
sys.path.append('../')
import numpy as np
from my_packages import functions_new
from ms_entropy import FlashEntropySearch,FlashEntropySearchCore
from tqdm import trange,tqdm
from FPSim2.io import create_db_file
from FPSim2 import FPSim2Engine

def parse_args():
    parse = argparse.ArgumentParser(description='Passsing mass spectral algorithms') # Create parameter object
    # parse.add_argument('-a','--algorithm',type=str, help='Similarity algorithm of tandem mass matching used for library search'
    #                          , default='entropy') # add parameter to object
    parse.add_argument('-i', '--istype', type=str,
                       help='type of isdb'
                       , default='is0')  # add parameter to object
    args = parse.parse_args() # parse parameter object to get parse object
    return args


def convert_numpy(obj):
    '''
    Recursively convert numpy types to native Python types
    Avoiding failed json_dump
    '''
    if isinstance(obj, dict):
        return {k: convert_numpy(v) for k, v in obj.items()}
    elif isinstance(obj, list):
        return [convert_numpy(i) for i in obj]
    elif isinstance(obj, np.ndarray):
        return obj.tolist()
    elif isinstance(obj, (np.integer, np.int64, np.uint32)):
        return int(obj)
    elif isinstance(obj, (np.floating, np.float32)):
        return float(obj)
    else:
        return obj


if __name__ == '__main__':
    t = time.time()

    args = parse_args()

    '''
    Load library file and create FP database
    '''

    library_map = {
        'is0': '0',
        'is1': '1',
        'is2': '2'
    }


    EL = library_map[args.istype]
    ID2SMILE = functions_new.json_load(f'../msdb/data/ICofIS/chem/is{EL}_Id2Smile.json')
    FPDB_FILE = f'../msdb/data/ICofIS/chem/is{EL}_fp_db.h5'
    OUTPUT_MATRIX = f'../msdb/data/ICofIS/chem/is{EL}_ICofIS.npy'
    IDX2SMILE = [[s['smile'], idx] for idx, s in enumerate(ID2SMILE)]

    FPDB = f'../msdb/data/ICofIS/chem/{args.istype}_fp_db.h5'
    if not os.path.exists(FPDB):
        create_db_file(
            mols_source=IDX2SMILE,
            filename=FPDB,
            mol_format='smiles',  # required
            fp_type='Morgan',
            fp_params={'radius': 2, 'fpSize': 1024}
        )


    '''
    Search 
    '''
    fpe = FPSim2Engine(fp_filename=FPDB)
    hit_dict = {}
    for idx in trange(len(ID2SMILE)):
        try:
            SMILE = ID2SMILE[idx]['smile']
            results = fpe.similarity(SMILE, threshold=0.75, metric='dice', n_workers=1)  # [(LIBidx,similarity1),(),...]
            results = [x for x in results if x[0] != idx] # Remove self comparison
            fpidx_list = [x[0] for x in results]
            fpsim_list = [x[1] for x in results]

            if fpidx_list:
                hit_dict[idx] = {'idx': fpidx_list, 'sim_list': fpsim_list}  # {LIBidx: {},{},..}
        except:print(idx)

    '''
    Save results
    '''
    ECHEM_MATRIX = np.full((len(IDX2SMILE), len(IDX2SMILE)), 0.0, dtype=float)
    for key, value in hit_dict.items():
        QueryIdx = key
        MatchIdxs = value['idx']
        MatchSim = value['sim_list']
        for idx, MatchIdx in enumerate(MatchIdxs):
            ECHEM_MATRIX[QueryIdx, MatchIdx] = MatchSim[idx]
    np.save(OUTPUT_MATRIX, ECHEM_MATRIX)

    print(f'Finished in {(time.time() - t) / 60:.2f} min')
