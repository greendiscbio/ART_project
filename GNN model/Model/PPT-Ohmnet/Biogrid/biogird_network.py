import requests
import json
import networkx as nx
import mygene
import os, io
import numpy as np
import pandas as pd

def get_biogrid(genes_file):

    request_url = "https://webservice.thebiogrid.org" + "/interactions"
    access_key = "27abe6af33db260960fc9af27947f29a"

    infile = open(genes_file, "r")
    geneList = infile.read().split("\n")
    infile.close()

    # These parameters can be modified to match any search criteria following
    # the rules outlined in the Wiki: https://wiki.thebiogrid.org/doku.php/biogridrest
    params = {
        "accesskey": access_key,
        "format": "json",  # Return results in TAB2 format
        "geneList": "XPC" ,#.join(geneList),  # Must be | separated
        "searchNames": "true",  # Search against official names
        'interSpeciesExcluded': 'true', # interactions w/ different species are excluded
        'throughputTag': 'high', 
        "includeInteractors": "true",  # set to false to get interactions between genes
        "taxId": 9606,  # Limit to Homo sapiens
        "includeInteractorInteractions": "false"
    }

    r = requests.get(request_url, params=params)
    interactions = r.json()
    print(interactions)
    # interactions.to_csv("GNN model/Model/PPT-Ohmnet/Biogrid/dataset.csv")

    # Create a hash of results by interaction identifier
    data = {}
    for interaction_id, interaction in interactions.items():
        data[interaction_id] = interaction
        # Add the interaction ID to the interaction record, so we can reference it easier
        data[interaction_id]["INTERACTION_ID"] = interaction_id

    # Load the data into a pandas dataframe
    dataset = pd.DataFrame.from_dict(data, orient="index")

    # Re-order the columns and select only the columns we want to see

    columns = [
        "INTERACTION_ID",
        "ENTREZ_GENE_A",
        "ENTREZ_GENE_B",
        "OFFICIAL_SYMBOL_A",
        "OFFICIAL_SYMBOL_B",
        "EXPERIMENTAL_SYSTEM",
        "PUBMED_ID",
        "PUBMED_AUTHOR",
        "THROUGHPUT",
        "QUALIFICATIONS",
    ]
    dataset = dataset[columns]
    dataset.to_csv("GNN model/Model/PPT-Ohmnet/Biogrid/dataset.csv")
    edgelist = dataset[['OFFICIAL_SYMBOL_A', 'OFFICIAL_SYMBOL_B']]
    edgelist = edgelist[edgelist['OFFICIAL_SYMBOL_A'] != edgelist['OFFICIAL_SYMBOL_B']] # remove self loops
    edgelist.to_csv('GNN model/Model/PPT-Ohmnet/Biogrid/AD_BioGrid_PPI.edgelist', sep='\t', index=False, header=None)

    G_frozen = nx.read_edgelist("GNN model/Model/PPT-Ohmnet/Biogrid/AD_BioGrid_PPI.edgelist")  
    G = nx.Graph(G_frozen)
    
    original = G.number_of_nodes()

    # Delete nodes from components with less than 5 nodes
    nodes_to_remove = []
    # for component in list(nx.connected_components(G)):
    #     if len(component)<5:
    #         for node in component:
    #             G.remove_node(node)

    largest = G.number_of_nodes()
    lost = original - largest
    lost_percent = round((lost/original), 4)

    print('Whole network:', original, 'nodes')
    print('Biggest connected component:', largest, 'nodes')
    print('Percentage of lost genes/nodes:', lost, f'({lost_percent*100}%)')
    return edgelist


a = get_biogrid('GNN model/Data/rna_30.txt')
print(a)
