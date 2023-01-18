import networkx as nx
import pandas as pd

initial_G = nx.read_edgelist("2023 Work/Graph networks/GNN+RN/Data/biogrid_minimum.edgelist")  
G = nx.Graph(initial_G)
print(G)

df = pd.read_csv("2023 Work/Graph networks/GNN+RN/Data/most_relevant_genes_biogrid_found.csv")
G_a = G.subgraph(df.columns)
print(G_a)

nx.wr.write_edgelist(G_a, "2023 Work/Graph networks/GNN+RN/Data/selected_nodes_graph_edgelist.edgelist")
initial_G = nx.read_edgelist("2023 Work/Graph networks/GNN+RN/Data/selected_nodes_graph_edgelist.edgelist")  
G = nx.Graph(initial_G)
print(G)