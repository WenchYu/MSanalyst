
'''

'''
import sys,re,os,math,time,heapq,spectral_entropy
sys.path.append('../')
import pandas as pd
import numpy as np
import networkx as nx
import spectrum_utils.spectrum as sus
from tqdm import tqdm
from collections import Counter
from my_packages import ms2tools,functions,functions_new
from multiprocessing import Process
from my_packages.peaktools import neutral_loss,modified_cosine

def refMN_generate(SPEC_THRESHOLD, PEAK_THRESHOLD, SPEC_SIM_MATRIX, IDXtoCCMSID, ALGORITHM, GNPS_INFO):
    '''

    :param SPEC_THRESHOLD:
    :param PEAK_THRESHOLD:
    :param SPEC_SIM_MATRIX:
    :param IDXtoCCMSID:
    :param ALGORITHM:
    :param GNPS_LIBRARY_FILE:
    :return:
    '''
    # gnps_info = functions_new.json_load(GNPS_LIBRARY_FILE)
    G = nx.MultiGraph()
    for idx1 in range(len(IDXtoCCMSID)):
        CCMSID = IDXtoCCMSID[f'{idx1}']
        G.add_node(CCMSID)

    if ALGORITHM == 'neutral_loss':
        N_PEAK_MATRIX = np.load(f'../msdb/data/SpecSimMatrix/neutral_loss_peak_matrix.npy')

    elif ALGORITHM == 'modified_cosine':
        N_PEAK_MATRIX = np.load(f'../msdb/data/SpecSimMatrix/neutral_loss_peak_matrix.npy')

    else:
        N_PEAK_MATRIX = np.load(f'../msdb/data/SpecSimMatrix/cosine_peak_matrix.npy')

    for idx1 in range(len(IDXtoCCMSID)):
        CCMSID1 = IDXtoCCMSID[str(idx1)]
        # spec, spectrum = functions_new.GNPS_info_format(GNPS_INFO, CCMSID1)

        for idx2 in range(idx1):
            CCMSID2 = IDXtoCCMSID[str(idx2)]
            n_peak = N_PEAK_MATRIX[idx1, idx2]
            spec_sim = SPEC_SIM_MATRIX[idx1, idx2]
            if spec_sim >= SPEC_THRESHOLD and n_peak >= PEAK_THRESHOLD:
                edge_attr = {f'{ALGORITHM}_similarity': spec_sim, 'shared_peaks': n_peak,
                             'edge_type': f'{ALGORITHM}'}
                G.add_edge(CCMSID1, CCMSID2, **edge_attr)


    MN_output = f'../msdb/data/SpecSimNetwork/MN_{ALGORITHM}_{SPEC_THRESHOLD}_{PEAK_THRESHOLD}.graphml'
    nx.write_graphml(G, MN_output)

def NetworkACC(graphml, matrix):
    '''

    :param graphml: graphml file
    :param matrix: similarity matrices
    :return:
    '''
    matrix_df = pd.read_csv(matrix,index_col=0)
    matrix_df.index = matrix_df.index.astype(str)  # convert index to str

    G = nx.read_graphml(graphml)
    CLUSTERs = [c for c in nx.connected_components(G) if len(c) > 1]
    CLUSTERs = sorted(CLUSTERs, key=lambda x: len(x), reverse=True)  # descending sorting
    singletons = [node for node in G.nodes() if G.degree(node) == 0]
    num_singleton = len(singletons)
    N_CLUSTERs = len(CLUSTERs)
    N_CLUSTER_NODEs = len(G.nodes) - num_singleton  # number of nodes in components

    components_acc = 0.0
    cluster_ave_acc = []
    for c in CLUSTERs:
        if len(c) > 1:
            subgraph = G.subgraph(c)
            edges = subgraph.edges()
            nodes = subgraph.nodes()
            n = len(nodes)  # num_nodes_k
            m = len(edges)  # number of edges
            edge_scores = 0.0
            for edge in edges:
                node1 = edge[0]
                # smile1 = info_df[info_df['name'] == node1].smiles.astype(str).values[0]
                # mol1 = chem.MolFromSmiles(smile1)
                node2 = edge[1]
                # smile2 = info_df[info_df['name'] == node2].smiles.astype(str).values[0]
                # mol2 = chem.MolFromSmiles(smile2)
                # score_edge = cheminfo_tools.MCS(mol1, mol2)
                score_edge = matrix_df.loc[node1,node2]
                if score_edge > 0:
                    edge_scores += score_edge
                else:
                    score_edge = matrix_df.loc[node2, node1]
                    edge_scores += score_edge

            component_acc = (edge_scores / m) * n
            cluster_ave_acc.append(edge_scores / m)
            components_acc += component_acc
    try:
        Network_acc = components_acc / N_CLUSTER_NODEs
    except:
        Network_acc = components_acc
    return Network_acc,cluster_ave_acc, N_CLUSTERs, N_CLUSTER_NODEs

def network_accuracy(GRAPHML_FILE, MATRIX_FILE, IDXtoCCMSID_FILE):
    '''

    :param GRAPHML_FILE:
    :param MATRIX_FILE:
    :param IDXtoCCMSID_FILE:
    :return:
    '''
    IDXtoCCMSID = functions_new.json_load(IDXtoCCMSID_FILE)
    CCMSIDtoIDX = {value: key for key, value in IDXtoCCMSID.items()}

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
            score_edge = MATRIX[int(CCMSIDtoIDX[node1]), int(CCMSIDtoIDX[node2])]

            edge_scores += score_edge

        component_acc = (edge_scores / m) * n
        cluster_ave_acc.append(edge_scores / m)
        COMPONENTS_ACC += component_acc
    # try:
    NETWORK_ACC = COMPONENTS_ACC / N_CLUSTER_NODEs
    # except: Network_acc = COMPONENTS_ACC
    return NETWORK_ACC, cluster_ave_acc, N_CLUSTERs, N_CLUSTER_NODEs

def N20(graphml):
    '''

    :param graphml: graphml file
    :return:
    '''
    G = nx.read_graphml(graphml)
    percentile20 = len(G.nodes) * 0.2
    clusters = [len(c) for c in nx.connected_components(G)]
    clusters_sorted = sorted(clusters, reverse=True)
    num_nodes=0
    N20 = 0
    for num in clusters_sorted:
        num_nodes +=num
        if num_nodes>percentile20:
            N20 = num
            break
    return N20

def n20(GRAPHML_FILE):
    '''

    :param graphml: graphml file
    :return:
    '''
    G = nx.read_graphml(GRAPHML_FILE)
    percentile20 = len(G.nodes) * 0.2
    clusters_sorted = sorted([len(c) for c in nx.connected_components(G)], reverse=True) # set key as len to sort by cluster size
    num_nodes=0
    N20 = 0
    for num in clusters_sorted:
        num_nodes +=num
        if num_nodes>percentile20:
            N20 = num
            break
    return N20

def RatioCCC(graphml_file, graphml_info='./data/std_info.csv'):
    info_df = pd.read_csv(graphml_info,dtype=str)
    G = nx.read_graphml(graphml_file)
    clusters = [c for c in nx.connected_components(G) if len(c) > 1]
    clusters= sorted(clusters, key=lambda x: len(x), reverse=True) # descending sorting
    num_clusters = len(clusters)
    correctly_classified_component = 0.0
    purities =[]
    for cluster in clusters:
        subgraph = G.subgraph(cluster)
        nodes = subgraph.nodes()
        n = len(nodes)  # number of molecules
        classes = []
        for node in nodes:
            try:
                superclass = info_df[info_df['name'] == node].superclass.astype(str).values[0]
                classes.append(superclass)
            except:
                print(type(node))
                print(info_df[info_df['name'] == node])
        counts = Counter(classes)
        most_common_element, max_count = counts.most_common(1)[0]
        purity = max_count / n
        purities.append(purity)
        if purity >= 0.7:
            correctly_classified_component += purity
    try:
        ratio_correctly_classified_component: float = correctly_classified_component / num_clusters
        return ratio_correctly_classified_component, purities
    except:
        return 0, purities

def DistanceFromMCS(x, y):
    '''
    distance from (14, 0.7995)
    :param x:
    :param y:
    :return:
    '''
    target_x = 14
    target_y = 0.7995

    distance = math.sqrt((x - target_x) ** 2 + (y - target_y) ** 2)
    return distance

def MNfiltering(args):
    '''
    curating the generated SpecSimNetwork by spectral similarity threshold and matched peaks
    :param graphml: graphml file
    :param threshold: spectral similarity threshold
    :param mps: number of matched peaks
    :return: Curated SpecSimNetwork
    '''

    G = nx.read_graphml(args.graphml)
    for u, v, data in list(G.edges(data=True)):
        if data['pair_similarity'] < 0.7 or data['matched_peaks'] < 5:
            G.remove_edge(u, v)

    output = re.sub(
        r"std_(.*?)_(\d+\.\d+)_(\d+)\.graphml",
    f"std_\\1_0.7_5.graphml",
        os.path.basename(args.graphml)
    )
    nx.write_graphml(G, args.output)

def RCCC(graphml,canopus_tsv):
    '''
    Ratio of correctly classified component (RCCC)
    :param graphml:
    :param canopus_tsv:
    :return:
    '''
    canopus_df = pd.read_csv(canopus_tsv, sep='\t')
    canopus_df['mappingFeatureId'] = canopus_df['mappingFeatureId'].astype(str)  # convert int64 to str
    # 'NPC#superclass','NPC#superclass Probability', 'NPC#class', 'NPC#class Probability','mappingFeatureId'

    G = nx.read_graphml(graphml)
    clusters = [c for c in nx.connected_components(G) if len(c) > 1]
    num_clusters = len(clusters)
    correctly_classified_component = 0.0
    for c in clusters:
        subgraph = G.subgraph(c)
        nodes = subgraph.nodes()
        n = len(nodes)  # number of molecules
        classes = []
        for node in nodes:
            try:
                superclass = canopus_df[canopus_df['mappingFeatureId'] == node]['NPC#superclass'].astype(str).values[0]
                classes.append(superclass)
            except: pass # Skip nodes without class prediction
                # print(node)
                # print(canopus_df[canopus_df['mappingFeatureId'] == node])
        counts = Counter(classes)
        most_common_element, max_count = counts.most_common(1)[0] # most frequent annotation type counts
        purity = max_count / n
        if purity >= 0.5:
            correctly_classified_component += 1
    try:
        ratio_correctly_classified_component: float = correctly_classified_component / num_clusters
        return ratio_correctly_classified_component
    except:
        return 0

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

def ratio_of_correctly_classified_component(GRAPHML_FILE, GNPS_LIBRARY_FILE):
    '''

    :param GRAPHML_FILE:
    :param GNPS_LIBRARY_FILE:
    :return:
    '''

    G = nx.read_graphml(GRAPHML_FILE)
    gnps_info = functions_new.json_load(GNPS_LIBRARY_FILE)

    CLUSTERs = [c for c in nx.connected_components(G) if len(c) > 1]
    N_CLUSTERs = len(CLUSTERs)
    CORRECTLY_CLASSIFIED_COMPONENT = 0.0

    for CLUSTER in CLUSTERs:
        subgraph = G.subgraph(CLUSTER)
        NODEIDs = subgraph.nodes()
        N_NODEs = len(NODEIDs)  # number of molecular features
        classes = []
        for NODEID in NODEIDs:
            superclass = gnps_info[NODEID]['np_classifier_superclass']
            classes.append(superclass)
        #         except: pass # Skip nodes without class prediction
        #             # print(node)
        #             # print(canopus_df[canopus_df['mappingFeatureId'] == node])
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

def CCCP(graphml,cluster,canopus_tsv):
    '''
    Correctly classified component purity (CCCP)
    :param graphml:
    :param cluster: for cluster in nx.connected_components(G)
    :param canopus_tsv:
    :return:
    '''
    canopus_df = pd.read_csv(canopus_tsv, sep='\t')  # load canopus result
    canopus_df['mappingFeatureId'] = canopus_df['mappingFeatureId'].astype(str)  # convert int64 to str
    # columns ['ClassyFire#superclass','NPC#superclass','NPC#superclass Probability', 'NPC#class', 'NPC#class Probability','mappingFeatureId']

    G = nx.read_graphml(graphml)  # load graphml
    subgraph = G.subgraph(cluster)  # load cluster in G
    nodes = subgraph.nodes()
    n = len(nodes)  # number of clustered molecular features
    classes = []
    for node in nodes:
        try:
            superclass = canopus_df[canopus_df['mappingFeatureId'] == node]['NPC#superclass'].astype(str).values[0]
            classes.append(superclass)
        except:
            pass  # Skip nodes without class prediction
        # print(node)
        # print(canopus_df[canopus_df['mappingFeatureId'] == node])
    counts = Counter(classes)
    most_common_element, max_count = counts.most_common(1)[0]
    purity = max_count / n
    return purity

def extract_cluster(GRAPHML):
    """
    graphml_file (str): path of graphml_file
    Returns:
            - clusters (idlist): [{node_id1,node_id2,node_id3},{},...]
            - cluster_nodes (idlist): [node_id1,node_id2,node_id3]
            - cluster_map (dict): {node_id: cluster_id}。
    """
    if '.graphml' in GRAPHML:
        G = nx.read_graphml(GRAPHML)
    else: G = GRAPHML
    clusters = [list(c) for c in nx.connected_components(G) if len(c) >= 1] # Cluster containing ≥ 2 nodes
    cluster_nodes = [node for node in G.nodes() if G.degree(node) > 0] # Cluster nodes
    cluster_map = {node: idx for idx, cluster in enumerate(clusters) for node in cluster} # Assign a specific identifier for each cluster {node_id: cluster_id}

    return clusters, cluster_nodes, cluster_map

def get_groups():
    '''
    get algorithm names by group
    :return:(tuple)
    '''
    group1 = ["weighted_dot_product","hellinger","ms_for_id","manhattan","absolute_value",
              "intersection","lorentzian","ruzicka","motyka","fidelity","bhattacharya_2",
              "matusita","bhattacharya_1","squared_chord","vicis_symmetric_chi_squared_3",
              "probabilistic_symmetric_chi_squared","harmonic_mean","unweighted_entropy",
              "improved_similarity","clark","spectral_contrast_angle","pearson_correlation",
              "dot_product","inner_product","whittaker_index_of_association","dot_product_reverse"]
    group2 = ["avg_l","entropy","roberts","baroni_urbani_buser","modified_cosine","neutral_loss",
              "jaccard","dice","peak_percentage","squared_euclidean","euclidean","penrose_shape"]
    group3 = ["chebyshev","ms_for_id_v1","canberra","divergence","wave_hedges","penrose_size",
              "mean_character","symmetric_chi_squared"]
    return group1, group2, group3

def get_optimal_threshold():


    return {'modified_cosine': 0.72, 'weighted_dot_product': 0.89,  'ms_for_id': 0.20, 'manhattan': 0.53,
            'absolute_value': 0.53, 'intersection': 0.53,  'ruzicka': 0.36, 'motyka': 0.53, 'fidelity': 0.72,
            'bhattacharya_2': 0.75, 'matusita': 0.475, 'bhattacharya_1': 0.765,  'vicis_symmetric_chi_squared_3': 0.62,
            'probabilistic_symmetric_chi_squared': 0.67, 'harmonic_mean': 0.67, 'unweighted_entropy': 0.70,
            'clark': 0.51, 'spectral_contrast_angle': 0.72, 'pearson_correlation': 0.845, 'cosine':0.72,'dot_product': 0.72,
            'inner_product': 0.1, 'dot_product_reverse': 0.99, 'avg_l': 0.175, 'entropy': 0.535, 'roberts': 0.285, 'ms_for_id_v1': 0.98,
            'baroni_urbani_buser': 0.73, 'neutral_loss': 0.72, 'jaccard': 0.54, 'dice': 0.52, 'peak_percentage': 0.71,
            'lorentzian': 0.52, 'wave_hedges': 0.055, 'whittaker_index_of_association': 0.11,'penrose_size': 0.23,'hellinger': 0.18}

            # {canberra': 0.05, 'divergence': 0.03,'chebyshev': 0.88,'penrose_shape': 0.84,'symmetric_chi_squared': 0.950,'squared_chord': 0.72,'euclidean': 0.84,'mean_character': 0.985
    #          'squared_euclidean': 0.975,'improved_similarity': 0.51}

    # return {'probabilistic_symmetric_chi_squared': 0.7022775000324369, 'lorentzian': 0.5319789524647249, 'dice': 0.5230705534881557,
    #         'entropy': 0.5397611856460571, 'harmonic_mean': 0.7022775018950821, 'vicis_symmetric_chi_squared_3': 0.6653143992366495,
    #         'ruzicka': 0.3709537070601039, 'jaccard': 0.354160826519861, 'manhattan': 0.541161535307765, 'intersection': 0.5411615371704102,
    #         'unweighted_entropy': 0.7205553913671418, 'motyka': 0.5411615361624182, 'roberts': 0.3057035145858097,
    #         'ms_for_id': 0.2984325817604636, 'absolute_value': 0.5411615370170715, 'avg_l': 0.22281826132287585,
    #         'pearson_correlation': 0.8508176359647193, 'dot_product': 0.7308166837281885, 'dot_product_reverse': 0.9919738182950542,
    #         'neutral_loss': 0.7509661912918091, 'baroni_urbani_buser': 0.6997271638597352, 'bhattacharya_2': 0.7650999783391313,
    #         'modified_cosine': 0.7261109566688538, 'chebyshev': 0.883170560002327, 'euclidean': 0.8499067209496554,
    #         'symmetric_chi_squared': 0.950330055440451, 'ms_for_id_v1': 0.9803177434644996, 'squared_euclidean': 0.9774720075839154,
    #         'penrose_shape': 0.8499067209496558, 'penrose_size': 0.23929592630080698, 'inner_product': 0.11959596104061809,
    #         'canberra': 0.06121517360855744, 'clark': 0.5276956232585585, 'hellinger': 0.18045191023432017,
    #         'divergence': 0.03282778186427615, 'fidelity': 0.7356368242349007, 'bhattacharya_1': 0.7755457176820686,
    #         'improved_similarity': 0.5276956232585585, 'matusita': 0.4858374015666401, 'mean_character': 0.9863081140210852,
    #         'spectral_contrast_angle': 0.7308166837281885, 'squared_chord': 0.7356368223722556, 'wave_hedges': 0.05392860545931022,
    #         'weighted_dot_product': 0.8998003743496245, 'whittaker_index_of_association': 0.1153577892558455, 'peak_percentage': 0.7132857142857143}

    # return {'probabilistic_symmetric_chi_squared': 0.7, 'lorentzian': 0.53, 'dice': 0.52, 'entropy': 0.54,
    #         'harmonic_mean': 0.7, 'vicis_symmetric_chi_squared_3': 0.67, 'ruzicka': 0.37, 'jaccard': 0.35,
    #         'manhattan': 0.54, 'intersection': 0.54, 'unweighted_entropy': 0.7, 'motyka': 0.54, 'roberts': 0.31,
    #         'ms_for_id': 0.3, 'absolute_value': 0.54, 'avg_l': 0.22, 'pearson_correlation': 0.7, 'dot_product': 0.7,
    #         'dot_product_reverse': 0.7, 'neutral_loss': 0.7, 'baroni_urbani_buser': 0.7, 'bhattacharya_2': 0.7,
    #         'modified_cosine': 0.7, 'chebyshev': 0.7, 'euclidean': 0.7, 'symmetric_chi_squared': 0.7, 'ms_for_id_v1': 0.7,
    #         'squared_euclidean': 0.7, 'penrose_shape': 0.7, 'penrose_size': 0.24, 'inner_product': 0.12, 'canberra': 0.06,
    #         'clark': 0.53, 'hellinger': 0.18, 'divergence': 0.03, 'fidelity': 0.7, 'bhattacharya_1': 0.7,
    #         'improved_similarity': 0.53, 'matusita': 0.49, 'mean_character': 0.7, 'spectral_contrast_angle': 0.7,
    #         'squared_chord': 0.7, 'wave_hedges': 0.05, 'weighted_dot_product': 0.7, 'whittaker_index_of_association': 0.12, 'peak_percentage': 0.7}

def NetworkEvaluation(graphml_files, matrix='./data/MCSmatrix.csv', classinfo = './data/std_info.csv'):
    '''

    :param graphml_files: a directory containing graphml files
    :param matrix: pairwise chemical similarity table
    :param classinfo: class info of molecular features
    :return: A .csv will be generated in the current folder with SpecSimNetwork topology evaluation result
    '''


    n20s = []
    naccs = []
    rccs = []
    names = []

    for graphml in tqdm(graphml_files, total=len(graphml_files)):
        try:
            names.append(os.path.basename(graphml).replace('.graphml', ''))
            acc = NetworkACC(graphml, matrix)[0]
            naccs.append(acc)
            n20 = N20(graphml)
            n20s.append(n20)
            rcc = RatioCCC(graphml, classinfo)[0]
            rccs.append(rcc)

        except:
            print(graphml)
    df = pd.DataFrame({'name': names, 'N20': n20s, 'NACC': naccs, 'RCCC': rccs})

    # output_name = names[0].split('_')[1]
    output_name = f'{os.path.dirname(graphml_files[0])}/evaluation_summary.csv'
    df.to_csv(output_name,index=None)

def parallel_clustering(algorithm, args):
    args.self_clustering_method = algorithm
    ms2tools.self_clustering(args)

def connectivity_screening(args):
    '''

    :param args:
    :return: graphml_files and merge_guide.csv
    '''
    result_folder = f'{args.output}/{os.path.splitext(os.path.basename(args.quant_file))[0]}_result'
    merge_folder = result_folder + '/merge/'
    annoated_graphml_file = os.path.join(result_folder,
                                         f'{os.path.splitext(os.path.basename(args.mgf_file))[0]}_{args.self_clustering_method}_{args.self_clustering_similarity}_{args.self_clustering_peaks}.graphml')
    benchmark_method = args.self_clustering_method

    '''Clustering by different algorithm'''
    algorithms = get_groups()[0] + get_groups()[1] + get_groups()[2]
    processes = []
    for algorithm in algorithms:
        p = Process(target=parallel_clustering, args=(algorithm, args))
        p.start()
        processes.append(p)
    for p in processes:
        p.join()

    '''Additional connectivity screening'''
    graphmls_for_merge = [merge_folder + g for g in os.listdir(merge_folder) if 'graphml' in g]

    benchmark_graphml_file = annoated_graphml_file
    benchmark_clusters, benchmark_cluster_nodes, benchmark_cluster_map = exCluster(benchmark_graphml_file)

    df_output = f'{result_folder}/{benchmark_method}_merge_guide.csv'  # Creating result dataframe
    if os.path.exists(df_output):
        df = pd.read_csv(df_output)
    else:
        df_column1 = [','.join(sorted(cluster)) for cluster in benchmark_clusters]  # Conver set to str
        df_column1.append('Additional clusters')
        df = pd.DataFrame({benchmark_method: df_column1})
        df.to_csv(df_output)

    df_value = []  # Additional linkage of molecular features recognized by different algorithms are recorded only once
    for filtered_graphml in graphmls_for_merge:
        filtered_method = os.path.basename(filtered_graphml).replace('.graphml', '')
        df[filtered_method] = np.nan

        intersection_cluster_idxs = set()  # Record the cluster of the filtered SpecSimNetwork (set)
        for benchmark_cluster_node in benchmark_cluster_nodes:
            filtered_clusters, filtered_cluster_nodes, filtered_cluster_map = exCluster(filtered_graphml)
            filtered_clusters_idxs = set(filtered_cluster_map.values())

            benchmark_cluster_idx = benchmark_cluster_map[benchmark_cluster_node]
            benchmark_cluster = benchmark_clusters[benchmark_cluster_idx]  # (set)
            benchmark_cluster_str = ','.join(sorted(benchmark_cluster))  # nodes in cluster (str)

            try:
                filtered_cluster_idx = filtered_cluster_map[benchmark_cluster_node]
                intersection_cluster_idxs.add(filtered_cluster_idx)
                filtered_cluster = filtered_clusters[filtered_cluster_idx]  # (set)
                if len(filtered_cluster) < 20:  # limit the size of cluster
                    non_intersection = filtered_cluster.difference(benchmark_cluster)
                    non_intersection_str = ','.join(sorted(non_intersection))  # nodes in cluster (str)

                    if non_intersection_str not in df_value:
                        if benchmark_cluster_str in df[benchmark_method].values:
                            idx = df[df[benchmark_method] == benchmark_cluster_str].index[0]
                            df.loc[idx, filtered_method] = non_intersection_str
                            df_value.append(non_intersection_str)
            except:
                pass  # print(benchmark_cluster_node) # Some nodes don't form cluster in other SpecSimNetwork

            other_clusters = []
            remaining_cluster_idxs = filtered_clusters_idxs - intersection_cluster_idxs
            for idx in remaining_cluster_idxs:
                remaining_cluster = filtered_clusters[idx]
                if len(remaining_cluster) < 20:
                    remaining_cluster_str = ','.join(sorted(remaining_cluster))
                    if remaining_cluster_str not in other_clusters:
                        other_clusters.append(remaining_cluster_str)
                        idx = df[df[benchmark_method] == 'Additional clusters'].index[0]
                        df.loc[idx, filtered_method] = f'{other_clusters}'.strip("[]")
    df = df.dropna(axis=1, how='all')  # Delete all columns with empty values
    df.to_csv(df_output, index=None)

def connectivity_filtering(graphml,threshold, mps):
    '''
    curating the generated SpecSimNetwork by spectral similarity threshold and matched peaks
    :param graphml: graphml file
    :param threshold: spectral similarity threshold
    :param mps: number of matched peaks
    :return: Curated SpecSimNetwork in the same directory with input SpecSimNetwork
    '''

    G = nx.read_graphml(graphml)
    for u, v, data in list(G.edges(data=True)):
        if data['pair_similarity'] <= threshold or data['matched_peaks'] < mps:
            G.remove_edge(u, v)

    output = re.sub(
        r"std_(.*?)_(\d+\.\d+)_(\d+)\.graphml",
    f"std_\\1_{threshold}_{mps}.graphml",
        os.path.basename(graphml)
    )
    nx.write_graphml(G, output)

def mn_curating(G, topk):
    '''
    Limit the size of the cluster by setting the number of neighbors (topk) allowed for a node

    :param G: Network graph
    :param topk: Maximum neighbors allowed for a node
    :return: An curated G
    '''
    node_ids = list(G.nodes())
    for node_id in node_ids:
        if len(G[node_id]) > topk:
            edges = list(G.edges(node_id)) # List all the edges
            # Keep the topK most similar neighbors of a node
            result = \
                heapq.nlargest(topk,
                               [(data.get('pair_similarity', 0), neighbor) for neighbor, data in G[node_id].items()])
            topk_edges = [t[1] for t in result]
            for edge in edges:
                if edge[1] not in topk_edges:
                    G.remove_edge(edge[0], edge[1])
    return G

def self_clustering(args):
    '''

    :param args: args.output, args.quant_file, args.mgf_file, args.self_clustering_similarity
    :return:
    '''
    parent_folder = f'{args.output}/{os.path.splitext(os.path.basename(args.quant_file))[0]}_result/'
    if not os.path.exists(parent_folder): # Make sure the parent folder exist
        os.makedirs(parent_folder)
    exp_info = functions.mgf_process(args.mgf_file) # Loading query '.mgf' file
    G = nx.MultiGraph()  # Creating undirected graph
    for i, (id1, pm1, charge1, spec1) in exp_info.iterrows():
        pm1 = float(pm1)
        node_attr = {'pepmass': pm1}
        G.add_node(id1, **node_attr)  # add nodes and attributes

    # Self clustering
    for i, (id1, pm1, charge1, spec1) in tqdm(exp_info.iterrows(), total=len(exp_info)):
        pm1 = float(pm1)
        charge1 = int(charge1)
        spec1 = spectral_entropy.clean_spectrum(spec1, max_mz=pm1 - 0.01, noise_removal=0.01)
        mz1 = np.array(spec1[:, 0], dtype=np.float64)
        intensity1 = np.array(spec1[:, 1], dtype=np.float64)
        spectrum1 = sus.MsmsSpectrum(identifier=id1, precursor_mz=pm1, precursor_charge=charge1
                                     , mz=mz1, intensity=intensity1).remove_precursor_peak(0.01, "Da")
        peaks1 = len(spec1)
        G.add_node(id1, **{'num_fragments': peaks1})
        for j, (id2, pm2, charge2, spec2) in exp_info.iloc[:i, ].iterrows():
            pm2 = float(pm2)
            charge2 = int(charge2)
            spec2 = spectral_entropy.clean_spectrum(spec2, max_mz=pm2 - 0.01, noise_removal=0.01)
            mz2 = np.array(spec2[:, 0], dtype=np.float64)
            intensity2 = np.array(spec2[:, 1], dtype=np.float64)
            spectrum2 = sus.MsmsSpectrum(identifier=id2, precursor_mz=pm2, precursor_charge=charge2, mz=mz2,
                                         intensity=intensity2).remove_precursor_peak(0.01, "Da")

            sim, mps, pp = 0.0, 0, 0.0
            if args.self_clustering_method == 'modified_cosine':
                try:
                    result = modified_cosine(spectrum1, spectrum2, fragment_mz_tolerance=0.02)
                    sim = result.score
                    mps = result.matches
                    pp = mps / peaks1
                except:
                    pass
            elif args.self_clustering_method == 'neutral_loss':
                try:
                    result = neutral_loss(spectrum1, spectrum2, fragment_mz_tolerance=0.02)
                    sim = result.score
                    mps = result.matches
                    pp = mps / peaks1
                except:
                    pass
            else:
                try:
                    sim = spectral_entropy.similarity(spec1, spec2, method=args.self_clustering_method, ms2_da=0.02)
                    result = modified_cosine(spectrum1, spectrum2, fragment_mz_tolerance=0.02)
                    mps = result.matches
                    pp = mps / peaks1
                except:
                    pass
            if sim >= args.self_clustering_similarity \
                    and mps >= args.self_clustering_peaks:
                edge_attr = {'pair_similarity': sim, 'matched_peaks': mps, 'peak_percentage': pp,
                             'edge_type': args.self_clustering_method}
                G.add_edge(id1, id2, **edge_attr)
    # Output
    G = mn_curating(G,args.top_k)
    print('Self clustering finished!')
    MN_file = os.path.join(parent_folder,
                           f'{os.path.splitext(os.path.basename(args.mgf_file))[0]}_{args.self_clustering_method}_{args.self_clustering_similarity}_{args.self_clustering_peaks}.graphml')
    nx.write_graphml(G, MN_file)


# Experimental
def load_tp_from_matrix(ALGO, OptSpecThreshold, SHARED_PEAKs):
    '''
    1. Load spectral and chemical similarity matrices
    2. top 1 scoring

    '''

    DIR = '../msdb/data/std/SearchMatrix/'
    ESPEC_MATRIX = np.load(os.path.join(DIR, f'E_{ALGO}.npy'))  # spectral sim SpecSimMatrix
    EMP_MATRIX = np.load(os.path.join(DIR, f'EMP.npy'))  # matched_ion SpecSimMatrix
    ECHEM_MATRIX = np.load(os.path.join(DIR, f'EChemSim.npy'))  # chemical sim SpecSimMatrix

    Res = []
    if SHARED_PEAKs:
        EMaxIdx = np.full(ESPEC_MATRIX.shape[0], 0, dtype=int)  # Lib id
        # QueryIdx = np.full(ESPEC_MATRIX.shape[0], 0, dtype=int) # Feature id
        EMaxSpecSim = np.full(ESPEC_MATRIX.shape[0], 0.0, dtype=float)
        EMaxChemSim = np.full(ESPEC_MATRIX.shape[0], 0.0, dtype=float)
        EMaxPeak = np.full(ESPEC_MATRIX.shape[0], 0.0, dtype=float)

        for idx, EMP_ARR in enumerate(EMP_MATRIX):
            cols = np.where(EMP_ARR >= 6)[0]  # <tuple> Min matched peaks >= 5
            if cols.any():  #
                FilColIdx = np.argmax(ESPEC_MATRIX[idx, cols])
                best_col = cols[FilColIdx]
                EMaxIdx[idx] = best_col
                EMaxSpecSim[idx] = ESPEC_MATRIX[idx, best_col]
                EMaxChemSim[idx] = ECHEM_MATRIX[idx, best_col]
                EMaxPeak[idx] = EMP_MATRIX[idx, best_col]

        ELables = [1 if x >= 0.35 else 0 for x in EMaxChemSim]  # Chem sim ≥ 0.75 is considered true positive
        Y_PRED = [1 if x >= OptSpecThreshold else 0 for x in EMaxSpecSim]  # Spectral sim ≥ optimal spectral threshold
        RES = []
        for idx in range(len(ELables)):
            if Y_PRED[idx] == 1 and ELables[idx] == 1:
                TP = 1
            else:
                TP = 0
            Res.append({'y_pred': Y_PRED[idx],
                        'label': ELables[idx],
                        'tp': TP,
                        'query_idx': idx,
                        'lib_idx': EMaxIdx[idx],
                        'spec_sim': EMaxSpecSim[idx],
                        'n_peak': EMaxPeak[idx],
                        'chem_sim': EMaxChemSim[idx],
                        })


    else:
        EMaxIdx = [np.argmax(arr) for arr in ESPEC_MATRIX]  # Top1 scoring
        EMaxSpecSim = [ESPEC_MATRIX[idx][EMaxIdx[idx]] for idx in
                       range(len(EMaxIdx))]  # Spectral similarity with top 1 match
        EMaxChemSim = [ECHEM_MATRIX[idx][EMaxIdx[idx]] for idx in
                       range(len(EMaxIdx))]  # Chemical dice simialrity with top 1 match
        EMaxPeak = [EMP_MATRIX[idx][EMaxIdx[idx]] for idx in range(len(EMaxIdx))]
        ELables = [1 if x >= 0.35 else 0 for x in EMaxChemSim]
        Y_PRED = [1 if x >= OptSpecThreshold else 0 for x in EMaxSpecSim]
        RES = []
        for idx in range(len(ELables)):
            if Y_PRED[idx] == 1 and ELables[idx] == 1:
                TP = 1
            else:
                TP = 0
            Res.append({'y_pred': Y_PRED[idx],
                        'label': ELables[idx],
                        'tp': TP,
                        'query_idx': idx,
                        'lib_idx': EMaxIdx[idx],
                        'spec_sim': EMaxSpecSim[idx],
                        'n_peak': EMaxPeak[idx],
                        'chem_sim': EMaxChemSim[idx],
                        })

    return Res

def load_tp_from_ISmatrix(ALGO, OptSpecThreshold, EL=0, SHARED_PEAKs=False):
    '''
    1. Load spectral and chemical similarity matrices
    2. top 1 scoring

    '''

    DIR = '../msdb/data/std/SearchMatrix/'
    ESPEC_MATRIX = np.load(os.path.join(DIR, f'IS{EL}_{ALGO}.npy'))  # spectral sim SpecSimMatrix
    EMP_MATRIX = np.load(os.path.join(DIR, f'IS{EL}_{ALGO}_MP.npy'))  # matched_ion SpecSimMatrix
    ECHEM_MATRIX = np.load(os.path.join('../msdb/data/std/chem/', f'IS{EL}_ChemSim.npy'))  # chemical sim SpecSimMatrix

    Res = []
    if SHARED_PEAKs:
        EMaxIdx = np.full(ESPEC_MATRIX.shape[0], 0, dtype=int)  # Lib id
        # QueryIdx = np.full(ESPEC_MATRIX.shape[0], 0, dtype=int) # Feature id
        EMaxSpecSim = np.full(ESPEC_MATRIX.shape[0], 0.0, dtype=float)
        EMaxChemSim = np.full(ESPEC_MATRIX.shape[0], 0.0, dtype=float)
        EMaxPeak = np.full(ESPEC_MATRIX.shape[0], 0.0, dtype=float)

        for idx, EMP_ARR in enumerate(EMP_MATRIX):
            cols = np.where(EMP_ARR >= 6)[0]  # <tuple> Min matched peaks >= 5
            if cols.any():  #
                FilColIdx = np.argmax(ESPEC_MATRIX[idx, cols])
                best_col = cols[FilColIdx]
                EMaxIdx[idx] = best_col
                EMaxSpecSim[idx] = ESPEC_MATRIX[idx, best_col]
                EMaxChemSim[idx] = ECHEM_MATRIX[idx, best_col]
                EMaxPeak[idx] = EMP_MATRIX[idx, best_col]

        ELables = [1 if x >= 0.75 else 0 for x in EMaxChemSim]  # Chem sim ≥ 0.75 is considered true positive
        Y_PRED = [1 if x >= OptSpecThreshold else 0 for x in EMaxSpecSim]  # Spectral sim ≥ optimal spectral threshold
        RES = []
        for idx in range(len(ELables)):
            if Y_PRED[idx] == 1 and ELables[idx] == 1:
                TP = 1
            else:
                TP = 0
            Res.append({'y_pred': Y_PRED[idx],
                        'label': ELables[idx],
                        'tp': TP,
                        'query_idx': idx,
                        'lib_idx': EMaxIdx[idx],
                        'spec_sim': EMaxSpecSim[idx],
                        'n_peak': EMaxPeak[idx],
                        'chem_sim': EMaxChemSim[idx],
                        })


    else:
        EMaxIdx = [np.argmax(arr) for arr in ESPEC_MATRIX]  # Top1 scoring
        EMaxSpecSim = [ESPEC_MATRIX[idx][EMaxIdx[idx]] for idx in
                       range(len(EMaxIdx))]  # Spectral similarity with top 1 match
        EMaxChemSim = [ECHEM_MATRIX[idx][EMaxIdx[idx]] for idx in
                       range(len(EMaxIdx))]  # Chemical dice simialrity with top 1 match
        EMaxPeak = [EMP_MATRIX[idx][EMaxIdx[idx]] for idx in range(len(EMaxIdx))]
        ELables = [1 if x >= 0.35 else 0 for x in EMaxChemSim]
        Y_PRED = [1 if x >= OptSpecThreshold else 0 for x in EMaxSpecSim]
        RES = []
        for idx in range(len(ELables)):
            if Y_PRED[idx] == 1 and ELables[idx] == 1:
                TP = 1
            else:
                TP = 0
            Res.append({'y_pred': Y_PRED[idx],
                        'label': ELables[idx],
                        'tp': TP,
                        'query_idx': idx,
                        'lib_idx': EMaxIdx[idx],
                        'spec_sim': EMaxSpecSim[idx],
                        'n_peak': EMaxPeak[idx],
                        'chem_sim': EMaxChemSim[idx],
                        })

    return Res

if __name__ == '__main__':
    t = time.time()

    print(f'Finished in {(time.time() - t) / 60:.2f} min')
