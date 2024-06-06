import networkx as nx

G_biogrid = nx.read_edgelist("2023 Work/Graph networks/Biogrid/Data/biogrid_minimum.edgelist")
print(G_biogrid)
genes_file = "2023 Work/Data/Preprocessed_data/Feature selection/SFS/LR_CV_0/Nuevos datos RNA/300_features.txt"
infile = open(genes_file, "r")
genes = infile.read().split("\n")
infile.close()
G = G_biogrid.subgraph(genes)
print(G)
print(G.nodes())