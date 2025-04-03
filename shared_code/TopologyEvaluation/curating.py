# -*- coding: utf-8 -*-
# @Time :2025/4/1 15:38
# @Auther :Yuwenchao
# @Software : PyCharm
'''
Network curating for
'''
import re
import os
import time
import networkx as nx

def MNCurating(graphml,threshold, mps):
    '''
    curating the generated network by spectral similarity threshold and matched peaks
    :param graphml: graphml file
    :param threshold: spectral similarity threshold
    :param mps: number of matched peaks
    :return: Curated network
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

if __name__ == '__main__':
    t = time.time()
    dir = './std_quant_result/0.1_1/'
    graphmls = [dir+g for g in os.listdir(dir)]


    graphml = dir+'std_clark_0.1_1.graphml'
    thresholds = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
    mpss = [1,2,3,4,5]
    # for graphml in graphmls:
    for threshold in thresholds:
        for mps in mpss:
            MNCurating(graphml,threshold,mps)


    print(f'Finished in {(time.time() - t) / 60:.2f} min')
