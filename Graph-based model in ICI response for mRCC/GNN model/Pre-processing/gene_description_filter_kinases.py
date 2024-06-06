import pandas as pd
import numpy as np

path ='GNN model\Data\mapping_results.xlsx'
data = pd.read_excel(path)
print('Total genes: '+str(len(data)))

# Genes with empty description
nulls = data.isnull().sum()
print('Genes with empty description: '+str(len(nulls)))
data = data[data['description'].notna()]
data = data.drop(data.loc[data['description'].str.contains('Uncharacterized protein')].index)
print('Genes with description and not uncharacterized: '+str(len(data)))

# Delete genes which  do not codify proteins
data = data.drop(data.loc[data['description'].str.contains('non-protein')].index)
print('Genes which codiy proteins: '+str(len(data)))

# Order by gene in order to drop repeated values
data = data.sort_values("external_gene_name")

# Delete duplicates
data = data.drop_duplicates(subset=['external_gene_name'])

# Resume
print('Final coding genes:')
print(len(data.external_gene_name.value_counts()))

# Save
data.to_excel('GNN model\Data\Coding_genes\Coding_genes.xlsx', index=False)

# Genes related to the cancer pathway
data2 = data.loc[data['description'].str.contains('cancer')]

# Resume
print('Final cancer pathway genes:')
print(len(data2.external_gene_name.value_counts()))

# Save
data2.to_excel('GNN model\Data\Cancer_proteins\Cancer_proteins.xlsx', index=False)

# Genes related with kinases
data = data.loc[data['description'].str.contains('kinase')]

# Resume
print('Final kinase genes:')
print(len(data.external_gene_name.value_counts()))

data.to_excel('GNN model\Data\Kinase_genes\Kinase_genes.xlsx', index=False)

# Genes related with protein kinases
data = data.loc[data['description'].str.contains('protein')]

# Resume
print('Final protein kinase genes:')
print(len(data.external_gene_name.value_counts()))

data.to_excel('GNN model\Data\Kinase_proteins\kinase_proteins.xlsx', index=False)
