
'''

'''
import sys
sys.path.append('./')
import os
import networkx as nx
from my_packages import config

# def mn_merging(args):
#     mn1_file = args.mn1_file
#     mn1_basename = os.path.basename(mn1_file).replace('.graphml', '')
#     mn1_G = nx.read_graphml(mn1_file)
#     node_ids = list(mn1_G.nodes())
#
#     mn2_file = args.mn2_file
#     mn2_basename = os.path.basename(mn2_file).replace('.graphml', '')
#     mn2_G = nx.read_graphml(mn2_file)
#
#     for src, dst, edge_attrs in mn1_G.edges(data=True):
#         try:
#             if not mn2_G.has_edge(src, dst) or (mn2_G[src][dst][0] != edge_attrs and mn2_G[src][dst][1] != edge_attrs):
#                 mn2_G.add_edge(src, dst, **edge_attrs)
#         except:
#             if not mn2_G.has_edge(src, dst) or mn2_G[src][dst][0] != edge_attrs:
#                 mn2_G.add_edge(src, dst, **edge_attrs)
#
#     nx.write_graphml(mn2_G, f'{args.output}/{mn1_basename}—{mn2_basename}.graphml')
def mn_merging(args):
    '''
    Merging networks constructed by different spectral algorithms
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
    args = config.args
    mn_merging(args)

