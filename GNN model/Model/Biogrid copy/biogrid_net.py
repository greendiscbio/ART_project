import requests
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt

def create_net(i, depth_level, geneList, g, found, initial_gene, branch):
    request_url = "https://webservice.thebiogrid.org" + "/interactions"
    access_key = "27abe6af33db260960fc9af27947f29a"

    # print('i: '+str(i))
    # print('Current gene: '+str(g))
    if i < depth_level and found == False: # recursive calls
        print('i: '+str(i))
        print('Current gene: '+str(g))
        print('Found?: '+str(found))
        
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
            print('Number of high interactions from '+str(g)+': '+str(len(interactions)))
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
            while e <len(edgelist) and found == False:
                if (edgelist['OFFICIAL_SYMBOL_A'][e].upper() != g.upper() and edgelist['OFFICIAL_SYMBOL_A'][e].upper() != initial_gene):
                    found, branch = create_net(i+1, depth_level, geneList, edgelist['OFFICIAL_SYMBOL_A'][e].upper(), found, initial_gene, branch)
                    # branch.append(edgelist.iloc[e]['OFFICIAL_SYMBOL_A'].upper())
                elif edgelist['OFFICIAL_SYMBOL_B'][e].upper() != g.upper() and edgelist['OFFICIAL_SYMBOL_B'][e].upper() != initial_gene:
                    found, branch = create_net(i+1, depth_level, geneList, edgelist['OFFICIAL_SYMBOL_B'][e].upper(), found, initial_gene, branch)
                    # branch.append(edgelist.iloc[e]['OFFICIAL_SYMBOL_B'].upper())

                if (edgelist.iloc[e]['OFFICIAL_SYMBOL_A'].upper() != g.upper() and edgelist.iloc[e]['OFFICIAL_SYMBOL_A'].upper() in geneList and edgelist.iloc[e]['OFFICIAL_SYMBOL_A'].upper() != initial_gene):
                    found = True
                    if len(branch) == 0 or branch[len(branch)-1]!=edgelist.iloc[e]['OFFICIAL_SYMBOL_A'].upper():
                        branch.append(edgelist.iloc[e]['OFFICIAL_SYMBOL_A'].upper())
                    branch.append(g)
                    return (True, branch)
                elif (edgelist.iloc[e]['OFFICIAL_SYMBOL_B'].upper() != g.upper() and edgelist.iloc[e]['OFFICIAL_SYMBOL_B'].upper() in geneList and edgelist.iloc[e]['OFFICIAL_SYMBOL_B'].upper() != initial_gene):
                    found = True
                    if len(branch) == 0 or branch[len(branch)-1]!=edgelist.iloc[e]['OFFICIAL_SYMBOL_B'].upper():
                        branch.append(edgelist.iloc[e]['OFFICIAL_SYMBOL_B'].upper())
                    branch.append(g)
                    return (True, branch)
                
                e = e+1
                if found == True:
                    branch.append(g)
                    print(branch)
        else:
            print("No connections found")
    return found, branch               

def get_biogrid(genes_file):

    infile = open(genes_file, "r")
    geneList = infile.read().split("\n")
    infile.close()
    print("Initial gene list")
    print(geneList)
    branches = []
    for g in geneList:
        found, branch = create_net(0, 2, geneList, g, False, g, [])
        # ------------------------------------------------------------------------------ #
        #   0: current depth level                                                       #
        #   2: max depth level                                                           #
        #   geneList: initial gene list where connections need to be found               #
        #   g: current gene of the geneList in which we are looking for connections      #
        #   True: do we have found a match of the initial geneList?                      #
        #   initial_gene: familiy gene (not include loops)                               #
        # ------------------------------------------------------------------------------ #
        branches.append(branch)
        print(branches)
        if found:
            print("Rama encontrada para "+ str(g)+' : '+str(branch))
        else:
            print("Rama no encontrada para "+str(g))

    with open("GNN model/Model/Biogrid copy/Data/branches.txt","a") as file:
        for branch in branches:
            file.write(str(branch)+"\n")
            for b in range (len(branch)-1):
                with open("GNN model/Model/Biogrid copy/Data/biogrid_minimum.edgelist","a") as file2:
                    file2.write(str(branch[b])+" "+str(branch[b+1])+"\n")
    print(branches)

    # Load network
    initial_G = nx.read_edgelist("GNN model/Model/Biogrid copy/Data/biogrid_minimum.edgelist")  
    G = nx.Graph(initial_G)
    # Needed to add this interacction since is the only interaction of GJA9 gene (low interaction)
    G.add_node("GJA9")
    G.add_edge("GJA9","GRB2")
    nodes_to_add = ["GPR155", "FPR2", "ARAP1", "RNPS1", "SRPK2", "SDR42E1"]
    G.add_nodes_from(nodes_to_add)
    G.add_edge("GPR155", "FPR2")
    G.add_edge("FPR2", "ARAP1")
    G.add_edge("ARAP1", "GRB2")
    G.add_edge("XPC", "RNPS1")
    G.add_edge("RNPS1", "SRPK2")
    G.add_edge("SRPK2", "SDR42E1")
    print(G)
    print()
    nx.write_edgelist(G, "GNN model/Model/Biogrid copy/Data/biogrid_minimum.edgelist")

    # Draw graph
    color_map = []
    for node in G:
        if str(node) in geneList and str(node) != "GRB2" or str(node) =="GJA9":
            color_map.append('orange')
        else: 
            color_map.append('blue')  
    mapping = {'USP31':'RP11-20G6.3'}
    G = nx.relabel_nodes(G, mapping)
    nx.draw(G, node_color=color_map, with_labels=True)
    plt.show()
a = get_biogrid('GNN model/Model/Biogrid copy/Data/rna_30.txt')

##################### RNA 30 included in Biogrid ########################
# DLGAP4
# GRB2
# GPR155
# IL25
# KLHL5
# LANCL1-AS1
# LEMD1
# PCMT1
# USP31
# SDR42E1
# TARBP2
# TRIM43B
# XPC
#########################################################################