
'''

'''
import sys
sys.path.append('./')
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
    print(gfiles)
    gfiles
    if len(gfiles) < 2:
        print("Need at least two molecular SpecSimNetwork")

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
            merged_G = mn_merging(merged_G, mni_G)  # Gradually merge with subsequent files
            graphml_name += '—' + mni_basename
        nx.write_graphml(merged_G, f'{args.output}/{graphml_name}.graphml')

def mn_merge(args):
    '''
    First, extract the DB nodes from both graphs (level == 'DB').
    Find "extra DBs" that appear only in mn2_G but not in mn1_G.
    For each extra DB:
    • If there is no node with the same name in mn1_G, move it directly to mn1_G and change its level to match its neighbor in mn2_G (Ematch or ISmatch);
    • If there is already a node with the same name in mn1_G, skip it.
    Also, copy the edges from mn2_G to Ematch/ISmatch.
    '''

    # basename = os.path.splitext(os.path.basename(args.spectra_file))[0]
    # parent_folder = f'{args.output}/'

    # args.mn1_file = os.path.join(parent_folder,f'{basename}_modified_cosine_modified_cosine_0.72_3.graphml')
    # args.mn1_file = os.path.join(parent_folder,f'MCS_0.7_3.graphml')
    # args.mn2_file = os.path.join(parent_folder,f'{basename}_entropy_modified_cosine_0.72_3.graphml')

    mn1_name = os.path.basename(args.mn1_file).replace('.graphml','')
    mn2_name = os.path.basename(args.mn2_file).replace('.graphml','')
    parent_folder = os.path.dirname(args.mn1_file) # output list

    # Load graphml
    mn1_G = nx.read_graphml(args.mn1_file)
    mn2_G = nx.read_graphml(args.mn2_file)

    # Locate extra db node
    EM1 = {n for n, d in mn1_G.nodes(data=True) if d.get('level') == 'Ematch'}
    EM2 = {n for n, d in mn2_G.nodes(data=True) if d.get('level') == 'Ematch'}
    extra_EM = EM2 - EM1  # 只在 mn2 里出现的 E match
    print(extra_EM)
    # Preprocess extra db node
    for em in extra_EM:
        for db in mn2_G.neighbors(em):
            if mn2_G.nodes[db].get('level') != 'DB':
                continue

            if db not in mn1_G:
                mn1_G.add_node(db, **mn2_G.nodes[db]) # 3-a 如果 mn1 中没有这个 DB 节点，先复制过去

                mn1_G.nodes[em].update(mn2_G.nodes[em])
                mn1_G.nodes[db].update(mn2_G.nodes[db]) # 3-b 把 DB 节点的属性更新成 mn2_G 里的最新值

                mn1_G.add_edge(db, em, **mn2_G[db][em]) # 3-c 把 mn2_G 里 (DB, Ematch) 这条边也搬到 mn1_G （若已存在则覆盖属性）


    ISM1 = {n for n, d in mn1_G.nodes(data=True) if d.get('level') == 'ISmatch'}
    ISM2 = {n for n, d in mn2_G.nodes(data=True) if d.get('level') == 'ISmatch'}
    extra_ISM = ISM2 - ISM1  # 只在 mn2 里出现的 E match
    for em in extra_ISM:
        for db in mn2_G.neighbors(em):
            if mn2_G.nodes[db].get('level') != 'DB':
                continue

            if db not in mn1_G:
                mn1_G.add_node(db, **mn2_G.nodes[db]) # 3-a 如果 mn1 中没有这个 DB 节点，先复制过去

                mn1_G.nodes[em].update(mn2_G.nodes[em])
                mn1_G.nodes[db].update(mn2_G.nodes[db]) # 3-b 把 DB 节点的属性更新成 mn2_G 里的最新值

                mn1_G.add_edge(db, em, **mn2_G[db][em]) # 3-c 把 mn2_G 里 (DB, Ematch) 这条边也搬到 mn1_G （若已存在则覆盖属性）
    # Save
    mn1_out_file = os.path.join(parent_folder,f'{mn1_name}——{mn2_name}')
    nx.write_graphml(mn1_G, mn1_out_file)

    # print(basename)
    # print('')

if __name__ == '__main__':
    t = time.time()
    args = config.args
    # args.mn1_file = './msdb/data/kutz/KutzOsmac_result/KutzOsmac_entropy_modified_cosine_0.72_3.graphml'
    # args.mn2_file =  './msdb/data/kutz/KutzOsmac_result/KutzOsmac_peak_percentage_modified_cosine_0.72_3.graphml'
    mn_merge(args)




