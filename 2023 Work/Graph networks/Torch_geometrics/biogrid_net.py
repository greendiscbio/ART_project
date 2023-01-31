import requests
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt

def create_net(i, depth_level, geneList, g, initial_gene, G):
    request_url = "https://webservice.thebiogrid.org" + "/interactions"
    access_key = "27abe6af33db260960fc9af27947f29a"

    if i < depth_level: # recursive calls
        # print('i: '+str(i))
        print('Current gene: '+str(g))
        
        ##################################### QUERY ####################################################
        params = {
            "accesskey": access_key,
            "format": "json",  # Return results in TAB2 format
            "geneList": g.upper() ,  # Must be | separated
            "searchNames": "true",  # Search against official names
            'interSpeciesExcluded': 'true', # interactions w/ different species are excluded
            'throughputTag': 'high', 
            "includeInteractors": "true",  # set to false to get interactions between genes
            "taxId": 9606, #10090,  # Limit to Homo sapiens
            "includeInteractorInteractions": "false"
        }

        r = requests.get(request_url, params=params)
        interactions = r.json()
        if r.status_code == 200 and len(interactions) > 0:
            print('Number of high interactions from '+str(g)+': '+str(len(interactions)) + " => " + str(initial_gene))
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
                    "OFFICIAL_SYMBOL_A",
                    "OFFICIAL_SYMBOL_B"
                ]
                dataset = dataset[columns]
                edgelist = dataset[['OFFICIAL_SYMBOL_A', 'OFFICIAL_SYMBOL_B']]
                edgelist = edgelist[edgelist['OFFICIAL_SYMBOL_A'] != edgelist['OFFICIAL_SYMBOL_B']] # remove self loops
                
            ############################################## QUERY END #######################################################
            # print('final edge list: '+str(edgelist))
            e = 0
            while e <len(edgelist):
                if (edgelist['OFFICIAL_SYMBOL_A'][e].upper() != g.upper() and edgelist['OFFICIAL_SYMBOL_A'][e].upper() != initial_gene):
                    G = create_net(i+1, depth_level, geneList, edgelist['OFFICIAL_SYMBOL_A'][e].upper(), initial_gene, G)
                    G.add_node(edgelist['OFFICIAL_SYMBOL_A'][e].upper())
                    G.add_edge(g, edgelist['OFFICIAL_SYMBOL_A'][e].upper())
                elif edgelist['OFFICIAL_SYMBOL_B'][e].upper() != g.upper() and edgelist['OFFICIAL_SYMBOL_B'][e].upper() != initial_gene:
                    G = create_net(i+1, depth_level, geneList, edgelist['OFFICIAL_SYMBOL_B'][e].upper(), initial_gene, G)
                    G.add_node(edgelist['OFFICIAL_SYMBOL_B'][e].upper())
                    G.add_edge(g, edgelist['OFFICIAL_SYMBOL_B'][e].upper())             
                e = e+1
        else:
            print("No connections found")
    return G               

def get_biogrid(genes_file):

    infile = open(genes_file, "r")
    geneList = infile.read().split("\n")
    infile.close()
    print("Initial gene list")
    print(geneList)
    G=nx.Graph()
    for g in geneList:
        print("Gene : "+ str(g))
        G = create_net(0, 2, geneList, g, g, G)
        # ------------------------------------------------------------------------------ #
        #   0: current depth level                                                       #
        #   2: max depth level                                                           #
        #   geneList: initial gene list where connections need to be found               #
        #   g: current gene of the geneList in which we are looking for connections      #
        #   initial_gene: familiy gene (not include loops)                               #
        #   G:  final Graph which contains every connection                              #
        # ------------------------------------------------------------------------------ #


    # Load network    
    # Needed to add this interacction since is the only interaction of GJA9 gene (low interaction)
    # G.add_node("GJA9")
    # G.add_edge("GJA9","GRB2")
    # print(G)
    # print()

    # Save graph
    nx.write_edgelist(G, "2023 Work/Graph networks/Torch_geometrics/Data/g.edgelist")

    # Draw graph
    # color_map = []
    # for node in G:
    #     if str(node) in geneList and str(node) != "GRB2" or str(node) =="GJA9":
    #         color_map.append('orange')
    #     else: 
    #         color_map.append('blue')  
    # mapping = {'USP31':'RP11-20G6.3'}
    # G = nx.relabel_nodes(G, mapping)
    # nx.draw(G, node_color=color_map, with_labels=True)
    # plt.show()
a = get_biogrid('2023 Work\Data/Preprocessed_data/Feature selection/SFS/LR_CV_0/Nuevos datos RNA/300_features.txt')
