
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
from my_packages import ms2tools,functions
from multiprocessing import Process
from my_packages.peaktools import neutral_loss,modified_cosine

def NetworkACC(graphml, matrix):
    '''

    :param graphml: graphml file
    :param matrix: similarity matrices
    :return:
    '''
    matrix_df = pd.read_csv(matrix,index_col=0)
    matrix_df.index = matrix_df.index.astype(str)  # convert index to str

    G = nx.read_graphml(graphml)
    clusters = [c for c in nx.connected_components(G) if len(c) > 1]
    clusters = sorted(clusters, key=lambda x: len(x), reverse=True)  # descending sorting
    singletons = [node for node in G.nodes() if G.degree(node) == 0]
    num_singleton = len(singletons)
    Total_num = len(G.nodes) - num_singleton  # number of nodes in components

    components_acc = 0.0
    cluster_ave_acc = []
    for c in clusters:
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
                # mol1 = Chem.MolFromSmiles(smile1)
                node2 = edge[1]
                # smile2 = info_df[info_df['name'] == node2].smiles.astype(str).values[0]
                # mol2 = Chem.MolFromSmiles(smile2)
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
        Network_acc = components_acc / Total_num
    except: Network_acc = components_acc
    return Network_acc,cluster_ave_acc

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
    curating the generated network by spectral similarity threshold and matched peaks
    :param graphml: graphml file
    :param threshold: spectral similarity threshold
    :param mps: number of matched peaks
    :return: Curated network
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

def CCCP(graphml,cluster,canopus_tsv):
    '''
    Correctly classified component purity (CCCP)
    :param graphml:
    :param cluster: nx.connected_components(G)
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

def exCluster(graphml_file):
    """
    graphml_file (str): path of graphml_file
    Returns:
            - clusters (list): [{node_id1,node_id2,node_id3},{},...]
            - cluster_nodes (list): [node_id1,node_id2,node_id3]
            - cluster_map (dict): {node_id: cluster_id}。
    """
    G = nx.read_graphml(graphml_file)
    clusters = [c for c in nx.connected_components(G) if len(c) > 1] # Cluster containing ≥ 2 nodes
    cluster_nodes = [node for node in G.nodes() if G.degree(node) > 0] # Cluster nodes
    cluster_map = {node: idx for idx, cluster in enumerate(clusters) for node in cluster} # Assign a specific identifier for each cluster {node_id: cluster_id}

    return clusters, cluster_nodes, cluster_map

def get_groups():
    '''
    get algorithm names by group
    :return:(tuple)
    '''
    group1 = ['probabilistic_symmetric_chi_squared','lorentzian','dice','entropy','harmonic_mean','vicis_symmetric_chi_squared_3','ruzicka','jaccard','manhattan','intersection','unweighted_entropy','motyka','roberts','ms_for_id','absolute_value','avg_l']
    group2 = ['pearson_correlation','dot_product','dot_product_reverse','neutral_loss','baroni_urbani_buser','bhattacharya_2','modified_cosine']
    group3 = ['chebyshev','euclidean','symmetric_chi_squared','ms_for_id_v1','squared_euclidean','penrose_shape']
    group4 = ['penrose_size','inner_product','canberra','clark','hellinger','divergence', 'fidelity_distance', 'hattacharya_1', 'improved_similarity_distance', 'matusita_distance', 'mean_character_distance', 'spectral_contrast_angle_distance', 'squared_chord_distance', 'wave_hedges_distance', 'weighted_dot_product_distance', 'whittaker_index_of_association_distance']
    return group1, group2, group3, group4

def NetworkEvaluation(graphml_files, matrix='./data/MCSmatrix.csv', classinfo = './data/std_info.csv'):
    '''

    :param graphml_files: a directory containing graphml files
    :param matrix: pairwise chemical similarity table
    :param classinfo: class info of molecular features
    :return: A .csv will be generated in the current folder with network topology evaluation result
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

        intersection_cluster_idxs = set()  # Record the cluster of the filtered network (set)
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
                pass  # print(benchmark_cluster_node) # Some nodes don't form cluster in other network

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
    curating the generated network by spectral similarity threshold and matched peaks
    :param graphml: graphml file
    :param threshold: spectral similarity threshold
    :param mps: number of matched peaks
    :return: Curated network in the same directory with input network
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


if __name__ == '__main__':
    t = time.time()


    print(f'Finished in {(time.time() - t) / 60:.2f} min')
