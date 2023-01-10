import pandas as pd
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt

def BFS_alg(graph, start, geneList, all_gene_list):
    explored = [] # Visted nodes
    queue = [[start]] 
     
    while queue:
        path = queue.pop(0)
        node = path[-1]
         
        # Check if the current node has not been visited yet and is included in RNA file
        if node not in explored and node in all_gene_list:
            neighbours = graph[node]
             
            # Iterate over the neighbours of the node
            for neighbour in neighbours:
                new_path = list(path)
                new_path.append(neighbour)
                queue.append(new_path)
                 
                # Check if the neighbour node is the goal: included in RNA list and not the same as the start node (self-loops)
                if neighbour in geneList and neighbour != start:
                    print("Shortest path = ", *new_path)
                    return new_path
            explored.append(node)
 
    # Nodes may be not connected
    print("Not connecting has been found")
    return

     
def BFS(RNA_file, RNA_30_file, String_file):

    # Load geneList which should be connected 
    infile = open(RNA_30_file, "r")
    geneList = infile.read().split("\n")
    infile.close()

    # Load all known genes
    infile = open(RNA_file, "r")
    all_genes = infile.read().split("\n")
    infile.close()

    # Load file where raw string connections are annotated
    infile = open(String_file, "r")
    results = infile.read().split("\n")
    infile.close()
    edge_list = []

    # Create edge list file
    edge_list_file = 'GNN model\Model\String\Data\string_edgelist.edgelist'
    outfile = open(edge_list_file, "a")
    for line in results:
        l = line.split("\t")
        edge_list.append([l[0],l[1]])
        outfile.write(l[0]+ '\t'+ l[1] + '\n')
    outfile.close()

    # Load network
    initial_G = nx.read_edgelist(edge_list_file)  
    G = nx.Graph(initial_G)
    print(G)
    print()

    # Graph using dictionaries
    d = nx.to_dict_of_dicts(G)

    # Create graph to save the solution
    minimum_graph = nx.Graph()

    for gene in geneList:
        if G.has_node(gene): # Some genes of the initial list may not be found by String. Consequently they are not included in the edge list
            path = BFS_alg(d, gene, geneList, all_genes)
            minimum_graph.add_nodes_from(path)
            # print(minimum_graph.number_of_nodes())
            
            for p in range (len(path)-1):
                if path[p] not in all_genes:
                    print("Gene " + path[p] + " is not included in RNA file")
                    print()
                minimum_graph.add_edge(path[p], path[p+1])
            # print(minimum_graph.number_of_edges())
        else:
            print("Gene " + gene + " is not include in the graph")

    m_g = nx.to_dict_of_dicts(minimum_graph)
    print("Minimun graph where genes are connected at least with other gene of the initial gene list:")
    print(minimum_graph)
    print()
    print(m_g)

    # Save graph
    nx.write_edgelist(minimum_graph, "GNN model/Model/String/Data/string_minimum_graph.edgelist")

    # Draw graph
    color_map = []
    for node in minimum_graph:
        if str(node) in geneList:
            color_map.append('orange')
        else: 
            color_map.append('blue')  
    nx.draw(minimum_graph, node_color=color_map, with_labels=True)
    plt.show()



########################################### REQUERIMENTS BEFORE EXECUTING #############################################
# Before you execute this python file you need to:
# 1. Go to https://string-db.org/
# 2. Using "Multiple protein" option insert the list of genes/proteins/transcripts you need to connect.
# 3. After selecting them, you should see an initial PPI network. If every element is not connected to at least other, 
#    you should add as nodes as needed to connect every element of the initial list using "More" option.
# 4. When your network is connected export it as tabular text output.
# 5. save it as string_output.tsv
########################################################################################################################

################ PARAMS TAKEN INTO ACCOUNT FOR FINDING CONNECTED NETWORK ##################
# Network type: full STRING network
# Active interaction sources: all
# Minimum required interaction score: 0.55
###########################################################################################

BFS('GNN model/Data/all_genes_db.txt', 'GNN model/Data/rna_30.txt', 'GNN model\Model\String\Data\string_output.tsv')