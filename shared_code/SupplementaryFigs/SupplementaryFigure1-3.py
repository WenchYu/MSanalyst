
'''
Generating netwokrs using in-house microbial NP library (./data/std.mgf and ./data/std_quant.mgf)
The output file 'evaluation_summary.csv' was used for Supplementary Figure 1-3 and supplementary table 3
'''
import os,sys, time
sys.path.append("../")
import networkx as nx
from my_packages import ms2tools,config,topology
from tqdm import tqdm,trange

def mn_merging(args):
    '''
    Merging networks
    :param args:
    :return:
    '''
    mn1_file = args.mn1_file
    mn1_basename = os.path.basename(mn1_file).replace('.graphml', '')
    mn1_G = nx.read_graphml(mn1_file)

    mn2_file = args.mn2_file
    mn2_basename = os.path.basename(mn2_file).replace('.graphml', '')
    mn2_G = nx.read_graphml(mn2_file)

    combined_G = nx.MultiGraph()
    combined_G.add_nodes_from(mn1_G.nodes(data=True))  # add nodes and edges from mn1_G
    combined_G.add_edges_from(mn1_G.edges(data=True))

    existing_edges = {}  # dick storing edges of mn1,{(node1,node2):[attibures1:xxx,...]}
    for u, v, data in mn1_G.edges(data=True):  # node1_id,node2_id,edge_attributes
        edge_key = (u, v)
        if edge_key not in existing_edges:
            existing_edges[edge_key] = []
        existing_edges[edge_key].append(data)

    for u, v, data in mn2_G.edges(data=True):
        edge_key = (u, v)
        if edge_key not in existing_edges:  # add edges not existing in mn1_G
            combined_G.add_edge(u, v, **data)
            existing_edges[edge_key] = [data]
        else:
            existing_data_list = existing_edges[edge_key]  # duplicate edges
            is_duplicate = False
            for existing_data in existing_data_list:
                if existing_data.items() == data.items():
                    is_duplicate = True
                    break
            if not is_duplicate:  # add edges with different attribtues
                combined_G.add_edge(u, v, **data)
                existing_edges[edge_key].append(data)
    nx.write_graphml(combined_G, f'{args.output}/{mn1_basename}—{mn2_basename}.graphml')  # output

if __name__ == '__main__':
    t = time.time()
    args = config.args
    args.quant_file = './data/std_quant.csv'
    args.mgf_file = './data/std.mgf'
    args.self_clustering_similarity = 0.7
    args.self_clustering_peaks = 5
    args.output = './networks/'
    if not os.path.exists(args.output):
        os.makedirs(args.output)

    '''Loading algorithms'''
    groups = topology.get_groups()
    group1_algorithms, group2_algorithms, group3_algorithms, group4_algorithms \
        = groups[0], groups[1], groups[2], groups[3]
    algorithms = groups[0] + groups[1] + groups[2] + groups[3]

    '''Generating networks'''
    for algorithm in tqdm(algorithms,total=len(algorithms)):
        args.self_clustering_method = algorithm
        ms2tools.self_clustering(args)

    '''Merging networks'''
    graphmls = [args.output + g for g in os.listdir(args.output)
                if any(algorithm in g for algorithm in group2_algorithms)] # group2_algorithms
    args.output = './0.7_5/group2_merging/' # './0.7_5/group2_merging/'
    if not os.path.exists(args.output):
        os.makedirs(args.output)
    for i in trange(len(graphmls)):
        args.mn1_file = graphmls[i]
        print(graphmls[i])
        for j in range(i):
            args.mn2_file = graphmls[j]
            mn_merging(args)

    '''Topology evaluation'''
    merged_graphmls = [args.output + g for g in os.listdir(args.output)]
    topology.NetworkEvaluation(merged_graphmls) # graphmls








    print(f'Finished in {(time.time() - t) / 60:.2f} min')
