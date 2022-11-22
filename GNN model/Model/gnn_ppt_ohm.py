import networkx as nx
import pandas as pd
import mygene

def get_snap(genes_file):

    G = nx.read_edgelist('GNN model\Data\Graph_ppt_ohm\PPT-Ohmnet_tissues-combined.edgelist', nodetype=int, data=(('tissue', str),))

    tissues_edgelist = pd.read_csv('GNN model/Data/Graph_ppt_ohm/PPT-Ohmnet_tissues-combined.edgelist', sep='\t')
    kidney_specific = tissues_edgelist[tissues_edgelist['tissue'] == 'kidney']

    kidney_specific.to_csv('GNN model/Data/Graph_ppt_ohm/PPT-Ohmnet_tissues-kidney.edgelist', sep='\t', index=False)
    G_kidney = nx.read_edgelist('GNN model/Data/Graph_ppt_ohm/PPT-Ohmnet_tissues-kidney.edgelist', nodetype=int, data=(('tissue', str),))

    # List of genes to search for
    infile = open(genes_file, "r")
    genes = infile.read().split("\n")
    infile.close()
    print(len(genes))

    # Genes in PPT-Ohmnet are Entrez IDs, it is necessary to convert them to gene Symbols.
    mg = mygene.MyGeneInfo()
    out = mg.querymany(genes, scopes='ensembl.gene', fields='entrezgene', species='human')
    print(out)
    entrezgenes = []
    mapping = {}
    for o in out:
        print(o)
        entrezgenes.append(int(o['entrezgene']))
        mapping[int(o['entrezgene'])] = o['query']

    A_kidney_frozen = G_kidney.subgraph(entrezgenes)
    A_kidney = nx.Graph(A_kidney_frozen)
    original = A_kidney.number_of_nodes()

    # Delete nodes from components with less than 5 nodes
    nodes_to_remove = []
    # for component in list(nx.connected_components(A_kidney)):
    #     if len(component)<5:
    #         for node in component:
    #             A_kidney.remove_node(node)

    # Remove self-loops
    A_kidney.remove_edges_from(list(nx.selfloop_edges(A_kidney)))

    largest = A_kidney.number_of_nodes()
    lost = original - largest
    lost_percent = round((lost/original), 4)

    print('Whole network:', original, 'nodes')
    print('Biggest connected component:', largest, 'nodes')
    print('Percentage of lost genes/nodes:', lost, f'({lost_percent*100}%)')

    A_kidney_relabeled = nx.relabel_nodes(A_kidney, mapping)
    nx.write_edgelist(A_kidney_relabeled, 'GNN model/Data/Graph_ppt_ohm/AD_SNAP_PPI_kidney.edgelist')

    return A_kidney_relabeled


get_snap('GNN model\Data\Graph_ppt_ohm\RNA_30.txt')