'''
This file creates a PPI network with the given depth level (level) of a given list of genes 
(genes_found.txt.txt) based on the Biogrid database and shows it in the .edgelist file 
(biogrid_found_genes_level_x.edgelist). It also creates one file with the n genes not known by Biogrid or 
which have no connection with any other gene (final_not_included_genes_level_x.txt).

As input parameters we need a list of genes to look for in Biogrid and as we found some connection problems 
during the long execution, we add a file called left_genes.txt which at first contains the same lines (genes) 
as genes_found.txt file. When the execution is interrupted due to some problem we update the left_file.txt 
by deleting the genes which have already been looked for. Even though we delete this genes they are taking 
into account as probable neighbors since they are included in the genes_found.txt file
'''

import requests
import pandas as pd
import networkx as nx

def create_net(i, depth_level, geneList, g, found, initial_gene, branch):
    request_url = "https://webservice.thebiogrid.org" + "/interactions"
    access_key = "27abe6af33db260960fc9af27947f29a"

    if i < depth_level and found == False: # recursive calls
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
            e = 0
            while e <len(edgelist) and found == False:
                if (edgelist['OFFICIAL_SYMBOL_A'][e].upper() != g.upper() and edgelist['OFFICIAL_SYMBOL_A'][e].upper() != initial_gene):
                    found, branch = create_net(i+1, depth_level, geneList, edgelist['OFFICIAL_SYMBOL_A'][e].upper(), found, initial_gene, branch)
                elif edgelist['OFFICIAL_SYMBOL_B'][e].upper() != g.upper() and edgelist['OFFICIAL_SYMBOL_B'][e].upper() != initial_gene:
                    found, branch = create_net(i+1, depth_level, geneList, edgelist['OFFICIAL_SYMBOL_B'][e].upper(), found, initial_gene, branch)
                if (edgelist.iloc[e]['OFFICIAL_SYMBOL_A'].upper() != g.upper() and edgelist.iloc[e]['OFFICIAL_SYMBOL_A'].upper() in geneList and edgelist.iloc[e]['OFFICIAL_SYMBOL_A'].upper() != initial_gene and g in geneList):
                    print("entro")
                    found = True
                    if len(branch) == 0 or branch[len(branch)-1]!=edgelist.iloc[e]['OFFICIAL_SYMBOL_A'].upper():
                        branch.append(edgelist.iloc[e]['OFFICIAL_SYMBOL_A'].upper())
                    branch.append(g)
                    print(g)
                    return (True, branch)
                elif (edgelist.iloc[e]['OFFICIAL_SYMBOL_B'].upper() != g.upper() and edgelist.iloc[e]['OFFICIAL_SYMBOL_B'].upper() in geneList and edgelist.iloc[e]['OFFICIAL_SYMBOL_B'].upper() != initial_gene and g in geneList):
                    found = True
                    print("entro2")
                    if len(branch) == 0 or branch[len(branch)-1]!=edgelist.iloc[e]['OFFICIAL_SYMBOL_B'].upper():
                        branch.append(edgelist.iloc[e]['OFFICIAL_SYMBOL_B'].upper())
                        print("ee")
                        print(edgelist.iloc[e]['OFFICIAL_SYMBOL_B'].upper())
                    branch.append(g)
                    print(g)
                    return (True, branch)
                
                e = e+1
                
                if found == True:
                    print(g)
                    branch.append(g)
    return found, branch               

def get_biogrid(genes_file, genes_left_file, level):

    infile = open(genes_file, "r")
    geneList = infile.read().split("\n")
    infile.close()
    print("Initial gene list")
    print(len(geneList))
    
    infile = open(genes_left_file, "r")
    genes_left = infile.read().split("\n")
    infile.close()
    print("Left gene list")
    print(len(genes_left))
    with open('2023 Work/Graph networks/Data/biogrid_found_genes_level_'+level+'.edgelist',"a") as file2:
        for g in genes_left:
            print("Initial gene: " + str(g))
            found, branch = create_net(0, int(level), geneList, g, False, g, [])
            # ------------------------------------------------------------------------------ #
            #   0: current depth level                                                       #
            #   level: max depth level                                                       #
            #   geneList: initial gene list where connections need to be found               #
            #   g: current gene of the geneList in which we are looking for connections      #
            #   True: do we have found a match of the initial geneList?                      #
            #   initial_gene: familiy gene (not include loops)                               #
            # ------------------------------------------------------------------------------ #
            if found:
                # print(branch)
                for b in range (len(branch)-1):
                    file2.write(str(branch[b])+" "+str(branch[b+1])+"\n")

    file2.close()
    initial_G = nx.read_edgelist('2023 Work/Graph networks/Data/biogrid_found_genes_level_'+level+'.edgelist')  
    G = nx.Graph(initial_G)
    with open('2023 Work/Graph networks/Data/final_not_included_genes_level_'+str(level)+'.txt',"a") as file3:
        for g in geneList:
            if g not in G.nodes():
               file3.write(str(g)+"\n")

level = -1
while int(level) < 1:
    level = input("How many levels do you want to look for? (level>=1) ")

a = get_biogrid('2023 Work/Graph networks/Data/genes_found.txt', '2023 Work/Graph networks/Data/left_genes.txt', level)

#-----------------------------------------------#
# Show principal characteristics of the graph   #
#-----------------------------------------------#
# initial_G = nx.read_edgelist('2023 Work/Graph networks/Data/HVGS/Biogrid/biogrid_found_genes_level_2_hvgs.edgelist')  
# # sub = nx.read_edgelist('2023 Work/Graph networks/Data/ANOVA/Biogrid/biogrid_found_genes_level_2_anova.edgelist')  
# G = nx.Graph(initial_G)
# # sub = nx.Graph(sub)
# print(G)
# # G = G.subgraph(sub.nodes())
# print(G)
# # print(nx.k_components(G))
# print(nx.number_connected_components(G))


#--------------------------------------------------------------------------------------------------#
# Check if any gene were not included in the initial list (genes whose sobreexpression is unknown) #
#--------------------------------------------------------------------------------------------------#

# infile = open("2023 Work/Graph networks/Data/genes_found.txt", "r")
# geneList = infile.read().split("\n")
# infile.close()
# print(len(geneList))
# print(len(G.nodes))
# for g in G.nodes:
#     # print(g)
#     if g not in geneList:
#         print(g)


#-----------------------------------------------------------------------------#
# Study wich genes where not included in the edgelist (no connecctions found) #
#-----------------------------------------------------------------------------#

# # Whole gene list
# infile = open("2023 Work/Graph networks/Data/HVGS/Biogrid/included_genes_hvgs.txt", "r")
# geneList = infile.read().split("\n")
# infile.close()

# # Achieved edgelist (genes which were connected)
# initial_G = nx.read_edgelist("2023 Work/Graph networks/Data/HVGS/Biogrid/biogrid_found_genes_level_2_hvgs.edgelist")
# G = nx.Graph(initial_G)

# # Save not connected genes
# with open("2023 Work/Graph networks/Data/ANOVA/Biogrid/not_included_genes_2_anova","a") as file3:
#     for g in geneList:
#         if g not in G.nodes():
#             file3.write(str(g)+"\n")
# file3.close()
