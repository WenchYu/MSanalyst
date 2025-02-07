# -*- coding: utf-8 -*-
# @Time :2023/5/4 21:44
# @Auther :Yuwenchao
# @Software : PyCharm
'''

'''
import os
import networkx as nx
from my_packages import config


if __name__ == '__main__':
    args = config.args
    mn1_file  = args.mn1_file
    mn1_basename = os.path.basename(mn1_file).replace('.graphml', '')
    mn1_G = nx.read_graphml(mn1_file)
    node_ids = list(mn1_G.nodes())

    mn2_file = args.mn2_file
    mn2_basename = os.path.basename(mn2_file).replace('.graphml', '')
    mn2_G = nx.read_graphml(mn2_file)

    for src, dst, edge_attrs in mn1_G.edges(data=True):
        try:
            if not mn2_G.has_edge(src, dst) or (mn2_G[src][dst][0] != edge_attrs and mn2_G[src][dst][1] != edge_attrs):
                mn2_G.add_edge(src, dst, **edge_attrs)
        except:
            if not mn2_G.has_edge(src, dst) or mn2_G[src][dst][0] != edge_attrs:
                mn2_G.add_edge(src, dst, **edge_attrs)

    nx.write_graphml(mn2_G, f'{args.output}/{mn1_basename}—{mn2_basename}.graphml')
