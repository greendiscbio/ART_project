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
    return found, branch               

def get_biogrid(genes_file):

    infile = open(genes_file, "r")
    geneList = infile.read().split("\n")
    infile.close()
    print("Initial gene list")
    print(geneList)
    with open("2023 Work/Graph networks/Spectral_clustering+CNN/Data/biogrid_foun_genes_level_2.edgelist","a") as file2:
        for g in geneList:
            print("Initial gene: " + str(g))
            found, branch = create_net(0, 3, geneList, g, False, g, [])
            # ------------------------------------------------------------------------------ #
            #   0: current depth level                                                       #
            #   2: max depth level                                                           #
            #   geneList: initial gene list where connections need to be found               #
            #   g: current gene of the geneList in which we are looking for connections      #
            #   True: do we have found a match of the initial geneList?                      #
            #   initial_gene: familiy gene (not include loops)                               #
            # ------------------------------------------------------------------------------ #
            if found:
                print(branch)
                for b in range (len(branch)-1):
                    file2.write(str(branch[b])+" "+str(branch[b+1])+"\n")

    file2.close()
    initial_G = nx.read_edgelist("2023 Work/Graph networks/Spectral_clustering+CNN/Data/biogrid_foun_genes_level_2.edgelist")  
    G = nx.Graph(initial_G)
    with open("Data/not_included_genes","a") as file3:
        for g in geneList:
            if g not in G.nodes():
               file3.write(str(g)+"\n")
    file3.close()
a = get_biogrid('2023 Work/Graph networks/Biogrid/Data/genes_found.txt')