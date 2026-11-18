import sys
import argparse
from Bio import SeqIO
import pandas as pd
import pickle
import igraph as ig
import leidenalg
import networkx as nx
import matplotlib.pyplot as plt

def get_core_graph_by_threshold(leiden_graph, threshold):
    leiden_core_graphs = []
    for graph in leiden_graph:
        vertices=graph['vertices']
        edges=graph['edges']
        G = nx.Graph()
        G.add_weighted_edges_from(edges)

        # Remove edges with weight below the threshold
        edges_to_remove = [(u, v) for u, v, w in G.edges(data=True) if w['weight'] < threshold]
        G.remove_edges_from(edges_to_remove)
        core_edges = list(G.edges())
        # Extract connected components
        connected_components = list(nx.connected_components(G))

        # Keep only the largest connected component
        largest_cc = max(connected_components, key=len)  # Get the largest component
        # Get all nodes
        core_vertices = []
        if len(list(largest_cc)) > 1:
            core_vertices = (largest_cc)
        # print("original",vertices )
        # print("core_vertices",core_vertices)
        # print("edges", edges)
        new_graph = {'vertices':vertices, 'edges':edges, 'core_vertices':core_vertices, 'core_edges':core_edges}
        leiden_core_graphs.append(new_graph)

    return leiden_core_graphs

parser = argparse.ArgumentParser()
parser.add_argument("-g", "--graph", required=True, help="leiden graph")
parser.add_argument("-o", "--output", required=True, default="output", help="output directory name")

def main():
    args = parser.parse_args()
    graph = args.graph
    output_path = args.output

    with open(graph, 'rb') as file:
        leiden_graph = pickle.load(file)
    threshold = 5
    leiden_core_graphs = get_core_graph_by_threshold(leiden_graph, threshold)

    filehandler = open(
        output_path + "/leiden_core_graphs_" + str(threshold) , "wb")
    pickle.dump(leiden_core_graphs, filehandler)
    print("all core graphs",len(leiden_core_graphs))
if __name__ == "__main__":
    main()