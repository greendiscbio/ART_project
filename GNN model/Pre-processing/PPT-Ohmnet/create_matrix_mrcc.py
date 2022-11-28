import pandas as pd
import numpy as np

num_genes = input('Number of genes for the matrix: ')
num_nodes = input('Number of nodes for the GNN: ')

path ='GNN model\Data\PPT-Ohmnet/Graph_PPT_Ohmn/AD_SNAP_PPI_kidney_'+str(num_genes)+'_genes_'+str(num_nodes)+'_nodes.edgelist'
data = pd.read_csv(path, delimiter=' ')
data = data[[data.columns[0], data.columns[1]]]
data = data.append(pd.Series([data.columns[0], data.columns[1]], index=[data.columns[0], data.columns[1]]), ignore_index=True)
print(data)
# print(data)
data.to_csv('GNN model/Data/PPT-Ohmnet/mRCC_big_pool/network_edges_mrcc_'+str(num_genes)+'_genes_'+str(num_nodes)+'_nodes.tsv', sep="\t")

path ='Data_preprocessing/Prediction PFS/RNA+Clinic joined/New/Clinical_data_and_RNA_total_Features_PFS.csv'
rna = pd.read_csv(path)
Y = [] # Target column
# For each entry I classified it by its PFS value.
for i in range (len(rna)):
    if rna.PFS[i]<3: # If PFS is lower than 3 months, I will consider it as NonResponder (NR)
        Y.append(0)
    else:
        Y.append(1)# If PFS is over 3 months, I will consider it as Responder (R)

data = data.sort_values(data.columns[0])
genes = np.concatenate((data[data.columns[0]], data[data.columns[1]]))
genes = np.unique(genes)

final_genes = pd.DataFrame(columns=genes)
not_in = []
for g in genes:
    print(g)
    if g in rna.columns:
        final_genes[g] = rna[g]
    else:
        not_in.append(g)

final_genes['Y']=Y
print(final_genes.head())
print('Genes not present in RNA Matrix: ' + str(len(not_in)))
print(not_in)
final_genes.to_csv('GNN model\Data/PPT-Ohmnet/mRCC_big_pool/mrcc_protein_matrix_'+str(num_genes)+'_genes_'+str(num_nodes)+'_nodes.csv')