# -*- coding: utf-8 -*-
# @Time :2025/7/28 00:31
# @Auther :Yuwenchao
# @Software : PyCharm
'''
Hqtof spectra match against In-Silico database
Filter experimental library search against the in-silico library
'''
import time,argparse,sys,random
sys.path.append('../')
import numpy as np
from my_packages import functions_new
from ms_entropy import FlashEntropySearch,FlashEntropySearchCore
from tqdm import trange,tqdm
from spectral_entropy import similarity

def parse_args():
    parse = argparse.ArgumentParser(description='Passsing mass spectral algorithms') # Create parameter object
    # parse.add_argument('-a','--algorithm',type=str, help='Similarity algorithm of tandem mass matching used for library search'
    #                          , default='entropy') # add parameter to object
    parse.add_argument('-i', '--istype', type=str,
                       help='type of isdb'
                       , default='is0')  # add parameter to object
    args = parse.parse_args() # parse parameter object to get parse object
    return args

if __name__ == '__main__':
    t = time.time()

    args = parse_args()

    if args.istype == 'is0':
        FS_IS_LIBRARY = functions_new.json_load('../msdb/FS_isdb_e0.json')
        is_search = FlashEntropySearch()
        FS_IS_LIBRARY = is_search.build_index(FS_IS_LIBRARY)

    if args.istype == 'is1':
        FS_IS_LIBRARY = functions_new.json_load('../msdb/FS_isdb_e1.json')
        is_search = FlashEntropySearch()
        FS_IS_LIBRARY = is_search.build_index(FS_IS_LIBRARY)

    if args.istype == 'is2':
        FS_IS_LIBRARY = functions_new.json_load('../msdb/FS_isdb_e2.json')
        is_search = FlashEntropySearch()
        FS_IS_LIBRARY = is_search.build_index(FS_IS_LIBRARY)


    # Load filtered library
    t = time.time()
    FS_E_LIBRARY_NR = functions_new.json_load('../msdb/FS_hqtof.json')
    e_nr_search = FlashEntropySearch()
    FS_E_LIBRARY_NR = e_nr_search.build_index(FS_E_LIBRARY_NR)


    # # 分成10份
    # N_SPEC = len(FS_E_LIBRARY_NR)
    # chunk_size = N_SPEC // 10
    # if N_SPEC % 10 != 0:
    #     chunk_size += 1  # 保证最后一份也包含进去

    # for part in range(10):
    #     start_idx = part * chunk_size
    #     end_idx = min((part + 1) * chunk_size, N_SPEC)
    #
    #     if start_idx >= N_SPEC:
    #         break  # 防止越界
        # for i in trange(start_idx, end_idx, desc=f'Processing part {part}'):


    SPEC_IS0_SIM_MATRIC = []
    SPEC_IS0_NPEAK_MATRIC = []
    SPEC_IS0_SIM_MATRIC1 = []
    for SPEC in tqdm(FS_E_LIBRARY_NR, total=len(FS_E_LIBRARY_NR)):
        PM = SPEC['precursor_mz']
        PEAKs = SPEC['peaks']

        SIM_ARRAY,PEAK_ARRAY = is_search.open_search(
            precursor_mz=PM,
            peaks=PEAKs,
            ms1_tolerance_in_da=0.02,
            ms2_tolerance_in_da=0.02,
            output_matched_peak_number=True)

        SPEC_IS0_SIM_MATRIC.append(SIM_ARRAY)
        SPEC_IS0_NPEAK_MATRIC.append(PEAK_ARRAY)

        # SIM_ARRAY1 = is_search.hybrid_search(
        #     precursor_mz=PM, peaks=PEAKs,
        #     ms1_tolerance_in_da=0.02, ms2_tolerance_in_da=0.02)
        # SIM_ARRAY1 = np.round(SIM_ARRAY1.astype(float), decimals=2)
        # SPEC_IS0_SIM_MATRIC1.append(SIM_ARRAY1)

        # 保存每部分
        np.save(f'../msdb/lib_matrix/entropy_{args.istype}_sim_emis.npy', SPEC_IS0_SIM_MATRIC)
        np.save(f'../msdb/lib_matrix/entropy_{args.istype}_peak_emis.npy', SPEC_IS0_NPEAK_MATRIC)
        # np.save(f'../msdb/lib_matrix/modified_cosine_{args.istype}_sim1_emis.npy', SPEC_IS0_SIM_MATRIC1)

    print("All parts saved.")
    print(f'Finished in {(time.time() - t) / 60:.2f} min')
