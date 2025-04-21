
'''

'''
import os,sys, time
sys.path.append("../")
import networkx as nx
from my_packages import ms2tools,config,topology
from tqdm import tqdm,trange

if __name__ == '__main__':
    t = time.time()
    args = config.args
    args.quant_file = './data/std_quant.csv'
    args.mgf_file = './data/std.mgf'
    args.self_clustering_similarity = 0.1
    args.self_clustering_peaks = 1
    args.output = './networks/'
    result_folder = f'{args.output}/{os.path.splitext(os.path.basename(args.quant_file))[0]}_result/'
    if not os.path.exists(result_folder):
        os.makedirs(result_folder)

    '''Loading algorithms'''
    groups = topology.get_groups()
    group4_algorithms = groups[3]

    '''Generating networks'''
    for algorithm in tqdm(group4_algorithms, total=len(group4_algorithms)):
        args.self_clustering_method = algorithm
        topology.self_clustering(args)

    '''Filtering spectral connectivity'''
    group4_graphmls = [result_folder+g for g in os.listdir(result_folder) if 'graphml' in g]
    thresholds = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
    mpss = [1, 2, 3, 4, 5]
    for graphml in group4_graphmls:
        for threshold in thresholds:
            for mps in mpss:
                topology.connectivity_filtering(graphml, threshold, mps)

    '''Topology evaluation'''
    group4_graphmls = [result_folder + g for g in os.listdir(result_folder) if 'graphml' in g]
    topology.NetworkEvaluation(group4_graphmls) # graphmls




    print(f'Finished in {(time.time() - t) / 60:.2f} min')
