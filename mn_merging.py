
'''

'''
import sys
sys.path.append('./')
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

    for src, dst, edge_attrs in mn1_G.edges(data=True):
        try:
            if not mn2_G.has_edge(src, dst) or (mn2_G[src][dst][0] != edge_attrs and mn2_G[src][dst][1] != edge_attrs):
                mn2_G.add_edge(src, dst, **edge_attrs)
        except:
            if not mn2_G.has_edge(src, dst) or mn2_G[src][dst][0] != edge_attrs:
                mn2_G.add_edge(src, dst, **edge_attrs)

    nx.write_graphml(mn2_G, f'{args.output}/{mn1_basename}—{mn2_basename}.graphml')

if __name__ == '__main__':
    args = config.args
    mn_merging(args)

    mn1_file = '/Users/hehe/Desktop/cbfe23/cbfe23_quant_result/cbfe23_neutral_loss_0.7_5.graphml'
    mn2_file = '/Users/hehe/Desktop/cbfe23/cbfe23_quant_result/cbfe23_fbmn.graphml'
    # mn1_file = args.mn1_file
    mn1_basename = os.path.basename(mn1_file).replace('.graphml', '')
    mn1_G = nx.read_graphml(mn1_file)
    node_ids = list(mn1_G.nodes())

    # mn2_file = args.mn2_file
    mn2_basename = os.path.basename(mn2_file).replace('.graphml', '')
    mn2_G = nx.read_graphml(mn2_file)

    # node_merging
    for node in set(mn1_G.nodes()) | set(mn2_G.nodes()):
        attrs1 = mn1_G.nodes[node] if node in mn1_G else {} # # 获取两个图中该节点的属性
        attrs2 = mn2_G.nodes[node] if node in mn2_G else {}
        # 合并属性，确保不重复
        combined_attrs = {**attrs1, **attrs2}
        # 添加节点及其属性到新图 G
        mn2_G.add_node(node, **combined_attrs)

    # edge_merging
    for src, dst, edge_attrs in mn1_G.edges(data=True):

        if not mn2_G.has_edge(src, dst): # check if edge in mn2 exist
            mn2_G.add_edge(src, dst, **edge_attrs)
        else:
            existing_edges = mn2_G.get_edge_data(src, dst) # get all attributes of this node
            # 假设 mn2_G 是 MultiGraph，需要检查所有可能的边属性
            # 对于普通图，existing_edges 是一个字典，可能只有一个键值对
            if isinstance(existing_edges, dict) and existing_edges:
                # 遍历所有可能的边属性
                for key, attrs in existing_edges.items():
                    if attrs == edge_attrs:
                        break
                else:
                    # 如果没有找到匹配的属性，则添加新边
                    mn2_G.add_edge(src, dst, **edge_attrs)
            else:
                # 如果 mn2_G 中的边属性结构不匹配，直接添加新边
                mn2_G.add_edge(src, dst, **edge_attrs)


    nx.write_graphml(G, f'{args.output}/{mn1_basename}—{mn2_basename}.graphml')
    # nx.write_graphml(mn2_G, f'/Users/hehe/Desktop/cbfe23/cbfe23_quant_result/{mn1_basename}—{mn2_basename}.graphml')
