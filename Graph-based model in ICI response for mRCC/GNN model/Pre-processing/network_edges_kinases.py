import pandas as pd
import numpy as np

path ='GNN model\Data\Kinase_proteins\string_interactions_kinase_proteins.tsv'
data = pd.read_csv(path, delimiter='\t')
data = data[['#node1', 'node2']]
data.to_csv('GNN model/Data/Kinase_proteins/network_edges_kinase_proteins.tsv', sep="\t")

path ='Data_preprocessing/Prediction PFS/RNA+Clinic joined/New/Clinical_data_and_RNA_total_Features_PFS.csv'
rna = pd.read_csv(path)
Y = [] # Target column
# For each entry I classified it by its PFS value.
for i in range (len(rna)):
    if rna.PFS[i]<3: # If PFS is lower than 3 months, I will consider it as NonResponder (NR)
        Y.append(0)
    else:
        Y.append(1)# If PFS is over 3 months, I will consider it as Responder (R)

data = data.sort_values("#node1")
genes = data['#node1'].unique()
final_genes = pd.DataFrame(columns=genes)
for g in genes:
    print(g)
    final_genes[g] = rna[g]

final_genes['Y']=Y
print(final_genes.head())

final_genes.to_csv("GNN model\Data\Kinase_proteins\Kinase_gene_matrix.csv")