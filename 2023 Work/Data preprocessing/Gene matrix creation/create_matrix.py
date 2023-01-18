import pandas as pd
import numpy as np

path = input('Insert file path .edgelist: ')
# path ='2023 Work/Graph networks/GNN+RN/Data/biogrid_minimum.edgelist'
data = pd.read_csv(path, delimiter=' ')
data = data[[data.columns[0], data.columns[1]]]
data = data.append(pd.Series([data.columns[0], data.columns[1]], index=[data.columns[0], data.columns[1]]), ignore_index=True)
print(data)

path ='2023 Work/Data/Preprocessed_data/clinic_and_RNA_data_raw.csv'
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
    next = False
    if g in rna.columns:
        final_genes[g] = rna[g]
    else:
        while not next:
            change = input('Gene "' + str(g)+ '" is not included in RNA file. Do yo want to modify its name? (y/n) ')
            if change == "y":
                new_gene = input('Enter new gene name: ')
                if new_gene in rna.columns:
                    final_genes[g] = rna[new_gene]
                    next = True
                else:
                    print('Gene "' + str(new_gene)+ '" is not included in RNA file.')
            elif change == "n":
                next = True
                not_in.append(g)

final_genes['Y']=Y
print(final_genes.head())
print('Genes not included in RNA file: ' + str(len(not_in)))
print(not_in)

name = input('Insert file name to save it (without extension): ')
final_genes.to_csv('2023 Work/Data/Preprocessed_data/Gene Matrix/' + str(name) + '.csv')