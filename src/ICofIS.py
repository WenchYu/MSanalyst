# -*- coding: utf-8 -*-
# @Time :2025/7/28 00:31
# @Auther :Yuwenchao
# @Software : PyCharm
'''
Inner consistency of In-Silico database
Calculate spectral pairs in ISDB from each other to determine their Inner consistency
'''
import json
import time,argparse,sys,random
import pandas as pd
sys.path.append('../')
import numpy as np
from my_packages import functions_new
from ms_entropy import FlashEntropySearch,FlashEntropySearchCore
from tqdm import trange,tqdm

def parse_args():
    parse = argparse.ArgumentParser(description='Passsing mass spectral algorithms') # Create parameter object
    # parse.add_argument('-a','--algorithm',type=str, help='Similarity algorithm of tandem mass matching used for library search'
    #                          , default='entropy') # add parameter to object
    parse.add_argument('-i', '--istype', type=str,
                       help='type of isdb'
                       , default='is0')  # add parameter to object
    args = parse.parse_args() # parse parameter object to get parse object
    return args

def random_select(LIST,seed=None):
    '''
    Random pick 1/10 elements from the list object
    :param LIST:
    :return:
    '''
    if seed is not None:
        random.seed(seed)
    subset_size = max(1, len(LIST) // 10) # select 1/10 of the original list element
    random_indices = random.sample(range(len(LIST)), subset_size)     # pick indices randomly
    FS_IS_LIBRARY_subset = [FS_IS_LIBRARY[i] for i in sorted(random_indices)]     # extract subsets

    return FS_IS_LIBRARY_subset


if __name__ == '__main__':
    t = time.time()

    args = parse_args()

    if args.istype == 'is0':
        FS_IS_LIBRARY = functions_new.json_load('../msdb/FS_isdb_e0.json')
        FS_IS_LIBRARY = random_select(FS_IS_LIBRARY, seed=33)
        is_search = FlashEntropySearch()
        FS_IS_LIBRARY = is_search.build_index(FS_IS_LIBRARY)

    if args.istype == 'is1':
        FS_IS_LIBRARY = functions_new.json_load('../msdb/FS_isdb_e1.json')
        FS_IS_LIBRARY = random_select(FS_IS_LIBRARY, seed=33)
        is_search = FlashEntropySearch()
        FS_IS_LIBRARY = is_search.build_index(FS_IS_LIBRARY)

    if args.istype == 'is2':
        FS_IS_LIBRARY = functions_new.json_load('../msdb/FS_isdb_e2.json')
        FS_IS_LIBRARY = random_select(FS_IS_LIBRARY, seed=33)
        is_search = FlashEntropySearch()
        FS_IS_LIBRARY = is_search.build_index(FS_IS_LIBRARY)


    print(f'Total number: {len(FS_IS_LIBRARY)}')

    '''Save relationship of randomly selected ISDB'''
    DATA = []
    for SPEC in FS_IS_LIBRARY:
        ID = SPEC['id']
        SMILE = SPEC['smile']
        DATA.append({'id': ID,'smile': SMILE})
    with open(f'../msdb/lib_matrix/{args.istype}_Id2Smile.json','w') as f:
        json.dump(DATA,f)



    SPEC_IS0_SIM_MATRIC = []
    SPEC_IS0_NPEAK_MATRIC = []
    SPEC_IS0_SIM_MATRIC1 = []
    for SPEC in tqdm(FS_IS_LIBRARY,total=len(FS_IS_LIBRARY)):
        PM = SPEC['precursor_mz']
        PEAKs = SPEC['peaks']

        SIM_ARRAY, NPEAK_ARRAY = is_search.open_search(
            precursor_mz=PM, peaks=PEAKs,
            ms1_tolerance_in_da=0.02, ms2_tolerance_in_da=0.02,
            output_matched_peak_number=True)
        SIM_ARRAY = np.round(SIM_ARRAY.astype(float), decimals=2)
        SPEC_IS0_SIM_MATRIC.append(SIM_ARRAY)
        SPEC_IS0_NPEAK_MATRIC.append(NPEAK_ARRAY)

        SIM_ARRAY1 = is_search.hybrid_search(
            precursor_mz=PM, peaks=PEAKs,
            ms1_tolerance_in_da=0.02, ms2_tolerance_in_da=0.02)
        SIM_ARRAY1 = np.round(SIM_ARRAY1.astype(float), decimals=2)
        SPEC_IS0_SIM_MATRIC1.append(SIM_ARRAY1)

    # 保存每部分
    np.save(f'../msdb/lib_matrix/entropy_{args.istype}_sim_icofis.npy', SPEC_IS0_SIM_MATRIC)
    np.save(f'../msdb/lib_matrix/entropy_{args.istype}_peak_icofis.npy', SPEC_IS0_NPEAK_MATRIC)
    np.save(f'../msdb/lib_matrix/modified_cosine_{args.istype}_sim1_icofis.npy', SPEC_IS0_SIM_MATRIC1)

    print("All parts saved.")
    print(f'Finished in {(time.time() - t) / 60:.2f} min')
