import networkx as nx
import pandas as pd
import numpy as np
import mygene

def get_snap(genes_file):
    G = nx.read_edgelist('GNN model\Data\PPT-Ohmnet\Graph_PPT_Ohmn/PPT-Ohmnet_tissues-combined.edgelist', nodetype=int, data=(('tissue', str),))

    tissues_edgelist = pd.read_csv('GNN model/Data/PPT-Ohmnet/Graph_PPT_Ohmn/PPT-Ohmnet_tissues-combined.edgelist', sep='\t')
    kidney_specific = tissues_edgelist[tissues_edgelist['tissue'] == 'kidney']

    kidney_specific.to_csv('GNN model/Data/PPT-Ohmnet/Graph_PPT_Ohmn/Second big pool/PPT-Ohmnet_tissues-kidney.edgelist', sep='\t', index=False)
    G_kidney = nx.read_edgelist('GNN model/Data/PPT-Ohmnet/Graph_PPT_Ohmn/Second big pool/PPT-Ohmnet_tissues-kidney.edgelist', nodetype=int, data=(('tissue', str),))
    rna_genes = np.loadtxt('GNN model/Data/all_genes_db.txt', dtype=str,  delimiter=',', skiprows=0)

    # List of genes to search for
    infile = open(genes_file, "r")
    genes = infile.read().split("\n")
    infile.close()
    print('Quering '+ str(len(genes)-1)+' genes.')

    # Genes in PPT-Ohmnet are Entrez IDs, it is necessary to convert them to gene Symbols.
    mg = mygene.MyGeneInfo()
    out = mg.querymany(genes, scopes='symbol', fields='entrezgene', species='human')
    print(out)
    entrezgenes = []
    no_entrezgenes = []
    mapping = {}
    for o in out:
        print(o)
        if 'entrezgene' in o and o['query'] in rna_genes: #['ACACA', 'ADRA2B', 'ATP5F1A', 'CCN1', 'CCN2', 'DACH1', 'H2AX', 'H4C1', 'H4C11', 'H4C12', 'H4C13', 'H4C14', 'H4C15', 'H4C2', 'H4C3', 'H4C4', 'H4C5', 'H4C6', 'H4C8', 'H4C9', 'HNF1B', 'MARCHF8', 'NSD2', 'ORAI1', 'PECAM1', 'PRKN', 'TWNK', 'TXNIP']:
            entrezgenes.append(int(o['entrezgene']))
            mapping[int(o['entrezgene'])] = o['query']
        else:
            no_entrezgenes.append(o['query'])
    print('Genes without entrezgene: ' + str(len(no_entrezgenes)))
    A_kidney_frozen = G_kidney.subgraph(entrezgenes)
    A_kidney = nx.Graph(A_kidney_frozen)
    original = A_kidney.number_of_nodes()

    # Delete nodes from components with less than 5 nodes
    nodes_to_remove = []
    for component in list(nx.connected_components(A_kidney)):
        if len(component)<2:
            for node in component:
                A_kidney.remove_node(node)

    # Remove self-loops
    A_kidney.remove_edges_from(list(nx.selfloop_edges(A_kidney)))

    largest = A_kidney.number_of_nodes()
    lost = original - largest
    lost_percent = round((lost/original), 4)

    print('Whole network:', original, 'nodes')
    print('Biggest connected component:', largest, 'nodes')
    print('Percentage of lost genes/nodes:', lost, f'({lost_percent*100}%)')
    A_kidney_relabeled = nx.relabel_nodes(A_kidney, mapping)
    nx.write_edgelist(A_kidney_relabeled, 'GNN model/Data/PPT-Ohmnet/Graph_PPT_Ohmn/Second big pool/AD_SNAP_PPI_kidney_'+str(len(genes)-1)+'_genes_'+str(largest)+'_nodes.edgelist')
    print(A_kidney_relabeled)
    return A_kidney_relabeled



path = 'GNN model\Data\PPT-Ohmnet/mRCC_big_pool/Second big pool/mrcc genes second big pool.xlsx'
data = pd.read_excel(path)
print('Initial number of genes: '+str(len(data)))
gdas = input('Minimum gene-disease association score: ')

data = data[data['Score_gda'] > float(gdas)]
genes = data.Gene
genes = genes.unique()
print(len(genes))
np.savetxt('GNN model/Data/PPT-Ohmnet/mRCC_big_pool/Second big pool/mrcc_gene_list_('+str(len(genes))+'_genes).txt', genes, delimiter=',', fmt='%s') 
get_snap('GNN model/Data/PPT-Ohmnet/mRCC_big_pool/Second big pool/mrcc_gene_list_('+str(len(genes))+'_genes).txt')