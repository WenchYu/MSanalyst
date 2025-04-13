
'''

'''
import os
import time
import argparse
import networkx as nx
from my_packages import ms2tools, config

def list_graphml_files(folder_path,keywords=None):
    '''
    return the absolute path of selected graphml files
    '''
    graphml_files = []
    for filename in os.listdir(folder_path):
        if filename.endswith(".graphml") and "—" not in filename:
            if keywords is not None:
                for keyword in keywords:
                    if filename in keyword:
                        graphml_files.append(os.path.join(folder_path, filename))
                        break
    return graphml_files

def mn_merging(mn1_G,mn2_G):
    # mn1_file = mn1_file
    # mn1_basename = os.path.basename(mn1_file).replace('.graphml', '')
    # mn1_G = nx.read_graphml(mn1_file)
    node_ids = list(mn1_G.nodes())

    # mn2_file = mn2_file
    # mn2_basename = os.path.basename(mn2_file).replace('.graphml', '')
    # mn2_G = nx.read_graphml(mn2_file)

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

    return mn2_G
    # nx.write_graphml(mn2_G, f'{args.output}/{mn1_basename}—{mn2_basename}.graphml')

def multimerging(args):
    gfiles = list_graphml_files(args.output,keywords=args.merge_list)

    if len(gfiles) < 2:
        print("Need at least two molecular networks")

    if len(gfiles) == 2:
        mn1_basename = os.path.basename(gfiles[0]).replace('.graphml', '')
        mn2_basename = os.path.basename(gfiles[1]).replace('.graphml', '')
        graphml_name = mn1_basename +'-'+ mn2_basename
        mn1_G = nx.read_graphml(gfiles[0])
        mn2_G = nx.read_graphml(gfiles[1])
        merged_G = mn_merging(mn1_G, mn2_G)
        nx.write_graphml(merged_G, f'{args.output}/{graphml_name}.graphml')

    else:
        mn1_basename = os.path.basename(gfiles[0]).replace('.graphml', '')
        mn2_basename = os.path.basename(gfiles[1]).replace('.graphml', '')
        mn1_G = nx.read_graphml(gfiles[0])
        mn2_G = nx.read_graphml(gfiles[1])
        merged_G = mn_merging(mn1_G, mn2_G)
        graphml_name = mn1_basename + '—' + mn2_basename
        for i in range(2, len(gfiles)):
            mni_G = nx.read_graphml(gfiles[i])
            mni_basename = os.path.basename(gfiles[i]).replace('.graphml', '')
            merged_G = mn_merging(merged_G, mni_G)  # 逐步与后续文件合并
            graphml_name += '—' + mni_basename
        nx.write_graphml(merged_G, f'{args.output}/{graphml_name}.graphml')

if __name__ == '__main__':
    t = time.time()
    args = config.args
    multimerging(args)

    # folder_path = '/Users/hehe/Desktop/example_quant_result/'
    # kwds=['/Users/hehe/Desktop/example_quant_result/example_neutral_loss_0.1_1.graphml'
    #     ,'/Users/hehe/Desktop/example_quant_result/example_neutral_loss_0.7_4.graphml']
    # args.output = folder_path
    # args.merge_list = kwds



