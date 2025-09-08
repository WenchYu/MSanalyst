# -*- coding: utf-8 -*-
# @Time :2025/8/5 06:25
# @Auther :Yuwenchao
# @Software : PyCharm
'''
Evaluation the in-silico database
1. confusion SpecSimMatrix of spectral pairs
2. network evaluation
usage: python ICofIS_evaluation.py -i is0 # or is1 or is2
'''
import os,time,argparse,sys,ujson
sys.path.append('../')
import numpy as np
import networkx as nx
from tqdm import trange,tqdm
from my_packages import functions_new,evaluation
from collections import Counter

def parse_args():
    parse = argparse.ArgumentParser(description='Passsing mass spectral algorithms') # Create parameter object
    # parse.add_argument('-a','--algorithm',type=str, help='Similarity algorithm of tandem mass matching used for library search'
    #                          , default='entropy') # add parameter to object
    parse.add_argument('-i', '--istype', type=str,
                       help='type of isdb'
                       , default='is0')  # add parameter to object
    args = parse.parse_args() # parse parameter object to get parse object
    return args

def rccc(GRAPHML_FILE, IDXtoFID_FILE, FS_SPECTRA):
    '''

    :param GRAPHML_FILE:
    :param GNPS_LIBRARY_FILE:
    :return:
    '''
    IDXtoFID = functions_new.json_load(IDXtoFID_FILE)
    FIDtoIDX = {FID: IDX for IDX, FID in IDXtoFID.items()}

    G = nx.read_graphml(GRAPHML_FILE)

    CLUSTERs = [c for c in nx.connected_components(G) if len(c) > 1]
    N_CLUSTERs = len(CLUSTERs)
    CORRECTLY_CLASSIFIED_COMPONENT = 0.0

    for CLUSTER in CLUSTERs:
        subgraph = G.subgraph(CLUSTER)
        NODEIDs = subgraph.nodes()
        N_NODEs = len(NODEIDs)  # number of molecular features
        classes = []
        for NODEID in NODEIDs:
            IDX = FIDtoIDX[NODEID]
            superclass = FS_SPECTRA[int(IDX)]['superclass'][0]
            classes.append(superclass)
        counts = Counter(classes)
        most_common_element, max_count = counts.most_common(1)[0]  # most frequent annotation type counts
        purity = max_count / N_NODEs
        if purity >= 0.7:
            CORRECTLY_CLASSIFIED_COMPONENT += 1
    try:
        RATION_CORRECTLY_CLASSIFIED_COMPONENT: float = CORRECTLY_CLASSIFIED_COMPONENT / N_CLUSTERs
    except:
        RATION_CORRECTLY_CLASSIFIED_COMPONENT = 0
    return RATION_CORRECTLY_CLASSIFIED_COMPONENT

def network_accuracy(GRAPHML_FILE, MATRIX_FILE, IDXtoID_FILE):
    '''

    :param GRAPHML_FILE:
    :param MATRIX_FILE:
    :param IDXtoCCMSID_FILE:
    :return:
    '''
    IDXtoID = functions_new.json_load(IDXtoID_FILE)
    IDtoIDX = {value: key for key, value in IDXtoID.items()}

    MATRIX = np.load(MATRIX_FILE)
    G = nx.read_graphml(GRAPHML_FILE)

    CLUSTERs = [c for c in nx.connected_components(G) if len(c) > 1]
    CLUSTERs = sorted(CLUSTERs, key=lambda x: len(x), reverse=True)  # descending sorting
    SINGLETONs = [node for node in G.nodes() if G.degree(node) == 0]
    N_SINGLETONs = len(SINGLETONs)
    N_CLUSTERs = len(CLUSTERs)
    N_CLUSTER_NODEs = len(G.nodes) - N_SINGLETONs  # number of nodes in components

    COMPONENTS_ACC, cluster_ave_acc = 0.0, []
    for CLUSTER in CLUSTERs:
        subgraph = G.subgraph(CLUSTER)
        nodes = subgraph.nodes()
        edges = subgraph.edges()
        n = len(nodes)  # num_nodes_k
        m = len(edges)  # number of edges
        edge_scores = 0.0
        for edge in edges:
            node1 = edge[0]
            node2 = edge[1]
            score_edge = MATRIX[int(IDtoIDX[node1]), int(IDtoIDX[node2])]

            edge_scores += score_edge

        component_acc = (edge_scores / m) * n
        cluster_ave_acc.append(edge_scores / m)
        COMPONENTS_ACC += component_acc
    try:
        NETWORK_ACC = COMPONENTS_ACC / N_CLUSTER_NODEs
    except:
        NETWORK_ACC = COMPONENTS_ACC
    return NETWORK_ACC, cluster_ave_acc, N_CLUSTERs, N_CLUSTER_NODEs

if __name__ == '__main__':
    t = time.time()

    args = parse_args()

    DIR = '../msdb/data/ICofIS/'
    '''Load resource files'''
    library_map = {'is0': '0','is1': '1','is2': '2'}
    EL = library_map[args.istype]
    SEARCH_ALGORITHMs = ['entropy']
    ID2SMILE_FILE = os.path.join(DIR, f'chem/is{EL}_Id2Smile.json')
    ID2SMILE = functions_new.json_load(ID2SMILE_FILE)
    IDX2SMILE = [[s['smile'], idx] for idx, s in enumerate(ID2SMILE)]

    CHEM_MATRIX_FILE = os.path.join(DIR,f'chem/is{EL}_ICofIS.npy')

    print(f'is{EL}_db')
    '''Network evaluation'''
    IDXtoID = {idx: s['id'] for idx, s in enumerate(ID2SMILE)}
    IDXtoID_FILE = os.path.join(DIR, 'IDXtoFID.json')
    if not os.path.exists(IDXtoID_FILE):
        with open(IDXtoID_FILE, 'w') as f:
            ujson.dump(IDXtoID, f)

    MN_DIR = os.path.join(DIR,'network/')
    MN_FILES = [os.path.join(MN_DIR, f) for f in os.listdir(MN_DIR) if 'graphml' in f]

    TOPO_RES = []
    for MN_F in MN_FILES:
        # try:
        mn = os.path.basename(MN_F).replace('.graphml', '')
        N20 = evaluation.n20(MN_F)
        NACC, cacc, nc, ncn = network_accuracy(MN_F, CHEM_MATRIX_FILE,IDXtoID_FILE)
        # NETWORK_ACC, cluster_ave_acc, N_CLUSTERs, N_CLUSTER_NODEs
        # print(mn,nc,ncn,NACC)
        print(f'The network is {mn}.')
        print(f'Network accuracy is {NACC}, {ncn} nodes form {nc} clusters.')
        print(f'The number of clusters (cluster_ave_acc >= 0.7) is {len([c for c in cacc if c >= 0.7])}.')


    '''top1 scoring'''
    ISM_RES = []
    ESPEC_MATRIX = np.load(f'../msdb/data/ICofIS/entropy_is{EL}_sim_icofis.npy')
    # EMP_MATRIX = np.load(f'../msdb/data/ICofIS/entropy_is{EL}_peak_icofis.npy')
    ECHEM_MATRIX = np.load(f'../msdb/data/ICofIS/chem/is{EL}_ICofIS.npy')

    OT = 0.54  # Optimal threshold

    # Initialize cumulative counts
    TNs, FPs, FNs, TPs = 0, 0, 0, 0
    non_zero_pairs = 0

    # Get the indices of the upper triangle (excluding the diagonal)
    triu_indices = np.triu_indices(ESPEC_MATRIX.shape[0], k=1)

    # Iterate through the unique pairs
    for i, j in trange(len(triu_indices[0])):
        row_idx = triu_indices[0][i]
        col_idx = triu_indices[1][j]

        # Get the single label and prediction for this pair
        e_label = 1 if ECHEM_MATRIX[row_idx, col_idx] >= 0.75 else 0
        y_prediction = 1 if ESPEC_MATRIX[row_idx, col_idx] >= OT else 0

        # Manually update confusion matrix counts
        if e_label == 0 and y_prediction == 0:
            TNs += 1
        elif e_label == 0 and y_prediction == 1:
            FPs += 1
        elif e_label == 1 and y_prediction == 0:
            FNs += 1
        elif e_label == 1 and y_prediction == 1:
            TPs += 1

        # Count non-zero pairs
        if ESPEC_MATRIX[row_idx, col_idx] > 0:
            non_zero_pairs += 1

    print(f'The shape of spectral SpecSimMatrix is {ESPEC_MATRIX.shape}.')
    print(f'The number of unique non-zero spectral pairs is {non_zero_pairs}.')
    print(f'In the unique spectral pairs, TN = {TNs}, FP = {FPs}, FN = {FNs}, TP = {TPs}.')

    print(f'Finished in {(time.time() - t) / 60:.2f} min')
