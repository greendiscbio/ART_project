import networkx as nx
import pandas as pd
import numpy as np
import mygene
import os
os.chdir('2023 Work/Graph networks/Graph2Vec/ppi_creation/ppt-ohmnet')

def get_snap(genes_file):
    G = nx.read_edgelist('Data/PPT-Ohmnet_tissues-combined.edgelist', nodetype=int, data=(('tissue', str),))
    print('Initial PPT-Ohmnet graph: ')
    print(G)
    tissues_edgelist = pd.read_csv('Data/PPT-Ohmnet_tissues-combined.edgelist', sep='\t')
    kidney_specific = tissues_edgelist[tissues_edgelist['tissue'] == 'kidney']
    
    # kidney_specific.to_csv('Data/PPT-Ohmnet_tissues-kidney.edgelist', sep='\t', index=False)
    G_kidney = nx.read_edgelist('Data/PPT-Ohmnet_tissues-kidney.edgelist', nodetype=int, data=(('tissue', str),))
    print('kidney graph: ')
    print(G_kidney)

    # List of genes to search for
    infile = open(genes_file, "r")
    genes = infile.read().split("\n")
    infile.close()
    print('Quering '+ str(len(genes))+' genes.')

    # Genes in PPT-Ohmnet are Entrez IDs, it is necessary to convert them to gene Symbols.
    mg = mygene.MyGeneInfo()
    out = mg.querymany(genes, scopes='symbol', fields='entrezgene', species='human')
    print(len(out))
    entrezgenes = []
    no_entrezgenes = []; que=[]
    mapping = {}
    for o in out:
        # print(o)
        # if o['query'] in ['IGHV1-12' ]:
        #     print(o)
        if 'entrezgene' in o:
            entrezgenes.append(int(o['entrezgene']))
            que.append(o['query'])
            mapping[int(o['entrezgene'])] = o['query']
        elif o['query'] not in que:
            # print(o)
            no_entrezgenes.append(o['query'])

    print('Genes without entrezgene: ' + str(len(no_entrezgenes)))
    print('Included: '+ str(len(entrezgenes)))
    A_kidney_frozen = G_kidney.subgraph(entrezgenes)
    print(A_kidney_frozen)
    # print(A_kidney_frozen.number_of_nodes())
    A_kidney = nx.Graph(A_kidney_frozen)
    original = A_kidney.number_of_nodes()
    print(A_kidney)
    # Delete nodes from components with less than 5 nodes
    exit()
    for component in list(nx.connected_components(A_kidney)):
        if len(component)<5:
            for node in component:
                A_kidney.remove_node(node)

    # Remove self-loops
    A_kidney.remove_edges_from(list(nx.selfloop_edges(A_kidney)))
    print(A_kidney)
    largest = A_kidney.number_of_nodes()
    lost = original - largest
    lost_percent = round((lost/original), 4)

    print('Whole network:', original, 'nodes')
    print('Biggest connected component:', largest, 'nodes')
    print('Percentage of lost genes/nodes:', lost, f'({lost_percent*100}%)')
    A_kidney_relabeled = nx.relabel_nodes(A_kidney, mapping)
    # nx.write_edgelist(A_kidney_relabeled, 'Data/AD_SNAP_PPI_kidney_total_RNA_nodes_Ave_Axi.edgelist')
    print(A_kidney_relabeled)
    return A_kidney_relabeled



get_snap('../Data/disgenet_found_genes_in_rna_total_Ave_Axi.txt')