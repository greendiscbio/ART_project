import requests
import pandas as pd


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
            if (edgelist['OFFICIAL_SYMBOL_A'][e].upper() != g.upper()):
                found, branch = create_net(i+1, depth_level, geneList, edgelist['OFFICIAL_SYMBOL_A'][e].upper(), found, initial_gene, branch)
            else:
                found, branch = create_net(i+1, depth_level, geneList, edgelist['OFFICIAL_SYMBOL_B'][e].upper(), found, initial_gene, branch)
            # print(e)
            # print(g)
            # print(found)
            # print('len' +str(len(branch)))
            # print(edgelist.iloc[e]['OFFICIAL_SYMBOL_A'].upper())
            if (edgelist.iloc[e]['OFFICIAL_SYMBOL_A'].upper() != g.upper() and edgelist.iloc[e]['OFFICIAL_SYMBOL_A'].upper() in geneList and edgelist.iloc[e]['OFFICIAL_SYMBOL_A'].upper() != initial_gene):
                found = True
                if len(branch) == 0 or branch[len(branch)-1]!=edgelist.iloc[e]['OFFICIAL_SYMBOL_A'].upper():
                    branch.append(edgelist.iloc[e]['OFFICIAL_SYMBOL_A'].upper())
                branch.append(g)
                # print("entro")
                # print(branch)
                return (True, branch)
            elif (edgelist.iloc[e]['OFFICIAL_SYMBOL_B'].upper() != g.upper() and edgelist.iloc[e]['OFFICIAL_SYMBOL_B'].upper() in geneList and edgelist.iloc[e]['OFFICIAL_SYMBOL_B'].upper() != initial_gene):
                found = True
                if len(branch) == 0 or branch[len(branch)-1]!=edgelist.iloc[e]['OFFICIAL_SYMBOL_B'].upper():
                    branch.append(edgelist.iloc[e]['OFFICIAL_SYMBOL_B'].upper())
                branch.append(g)
                # print("entro2")
                # print(branch)
                return (True, branch)
            
            e = e+1
            if found == True:
                # print("gen")
                # print(g)
                branch.append(g)
                # print(branch)
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
        #   1: max depth level                                                           #
        #   geneList: initial gene list where connections need to be found               #
        #   g: current gene of the geneList in which we are looking for connections      #
        #   True: do we have found a match of the initial geneList?                      #
        #   initial_gene: familiy gene (not include loops)                               #
        # ------------------------------------------------------------------------------ #
        branches.append(branch)
        if found:
            print("Rama encontrada para "+ str(g)+' : '+str(branch))
        else:
            print("Rama no encontrada para "+str(g))

    with open("GNN model/Model/PPT-Ohmnet/Biogrid/Data/branches.txt","a") as file:
        for branch in branches:
            file.write(str(branch)+"\n")
            for b in range (len(branch)-1):
                with open("GNN model/Model/PPT-Ohmnet/Biogrid/Data/biogrid_minimum.edgelist","a") as file2:
                    file2.write(str(branch[b])+" "+str(branch[b+1])+"\n")
    print(branches)
a = get_biogrid('GNN model/Data/rna_30.txt')
print(a)
