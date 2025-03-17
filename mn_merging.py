# -*- coding: utf-8 -*-
# @Time :2023/5/4 21:44
# @Auther :Yuwenchao
# @Software : PyCharm
'''

'''
import os
import networkx as nx
from my_packages import config

def mn_merging(args):
    mn1_file = args.mn1_file
    mn1_basename = os.path.basename(mn1_file).replace('.graphml', '')
    mn1_G = nx.read_graphml(mn1_file)
    node_ids = list(mn1_G.nodes())

    mn2_file = args.mn2_file
    mn2_basename = os.path.basename(mn2_file).replace('.graphml', '')
    mn2_G = nx.read_graphml(mn2_file)

    # node_merging
    for node in set(mn1_G.nodes()) | set(mn2_G.nodes()):
        attrs1 = mn1_G.nodes[node] if node in mn1_G else {}  # get node attributes from two graph
        attrs2 = mn2_G.nodes[node] if node in mn2_G else {}

        combined_attrs = {**attrs1, **attrs2}  # merge attibute to make sure no duplicate
        mn2_G.add_node(node, **combined_attrs)  # add node and attibutes to new graph

    # edge_merging
    for src, dst, edge_attrs in mn1_G.edges(data=True):

        if not mn2_G.has_edge(src, dst):  # check if edge in mn2 exist
            mn2_G.add_edge(src, dst, **edge_attrs)
        else:
            existing_edges = mn2_G.get_edge_data(src, dst)  # get all attributes of this node
            # Assuming mn2_G is a MultiGraph, all possible edge attributes need to be checked
            # For a normal graph, existing_edges is a dictionary, which may have only one key-value pair
            if isinstance(existing_edges, dict) and existing_edges:
                for key, attrs in existing_edges.items():
                    if attrs == edge_attrs:
                        break
                else:
                    mn2_G.add_edge(src, dst, **edge_attrs)  # if no edge attributes, find new attributes
            else:
                mn2_G.add_edge(src, dst, **edge_attrs)  # if edge attributes in mn2_G doesn't match，add new edges

    nx.write_graphml(mn2_G, f'{args.output}/{mn1_basename}—{mn2_basename}.graphml')

if __name__ == '__main__':
    args = config.args
    mn_merging(args)


