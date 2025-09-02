# -*- coding: utf-8 -*-
# @Time :2025/7/6 23:25
# @Auther :Yuwenchao
# @Software : PyCharm
'''
Clustering by spectral sim SpecSimMatrix
'''
import time,argparse,sys,os
sys.path.append('../')
import networkx as nx
import numpy as np
from matchms.importing import load_from_mgf
from ms_entropy.file_io import spec_file
from my_packages import functions_new

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

    DIR = '../msdb/data/ICofIS/'
    '''Load resource files'''
    library_map = {'is0': '0','is1': '1','is2': '2'}
    EL = library_map[args.istype]
    # FS_IS_LIBRARY = functions_new.json_load(f'../msdb/FS_isdb_e{EL}.json')

    # 1/10 selected ISDB, [{id:,'smile':},{},{},...]
    ID2SMILE = functions_new.json_load(os.path.join(DIR,f'/chem/is{EL}_Id2Smile.json'))
    IDX2SMILE = [[s['smile'], idx] for idx, s in enumerate(ID2SMILE)]
    SPEC_SIM_MATRIX = np.load(os.path.join(DIR,f'/entropy_is{EL}_sim_icofis.npy'))


    '''Networking'''
    SPEC_THRESHOLD = 0.54 # optimal entropy threshold
    G = nx.MultiGraph()
    for idx, FS_SPECTRUM in enumerate(ID2SMILE):
        FID = FS_SPECTRUM['id']
        # PM = FS_SPECTRUM['precursor_mz']
        SMILE = FS_SPECTRUM['smile']
        NODE_ATTR = {'feature_id': FID, 'smile': SMILE}
        G.add_node(FID, **NODE_ATTR)

    for idx1 in range(len(ID2SMILE)):
        FID1 = ID2SMILE[idx1]['id']
        for idx2 in range(idx1):
            FID2 = ID2SMILE[idx2]['id']
            SPEC_SIM = SPEC_SIM_MATRIX[idx1, idx2]
            if SPEC_SIM >= SPEC_THRESHOLD:  # and NPEAK >= NPEAK_THRESHOLD:
                EDGE_ATTR = {'spectral_similairty': SPEC_SIM}
                G.add_edge(FID1, FID2, **EDGE_ATTR)
    nx.write_graphml(G, os.path.join(DIR,f'/network/entropy_{SPEC_THRESHOLD}.graphml'))
    print(f'Finished in {(time.time() - t) / 60:.2f} min')
