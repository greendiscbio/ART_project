from unittest import case
import pandas as pd
import numpy as np

path ='GNN model\Data\mapping_results.xlsx'
data = pd.read_excel(path)
print('Total genes: '+str(len(data)))

# Genes with empty description
data = data[data['description'].notna()]
data = data.drop(data.loc[data['description'].str.contains('Uncharacterized protein')].index)
print('Genes with description and not uncharacterized: '+str(len(data)))


# Order by gene in order to drop repeated values
data = data.sort_values("external_gene_name")

# Delete duplicates
data = data.drop_duplicates(subset=['external_gene_name'])

# Genes related to T cells
data2 = data.loc[data['description'].str.contains('T-cell')]

# Resume
print('Final T cell:')
print(len(data2.external_gene_name.value_counts()))


# Genes related with programmed cell death 1 
data3 = data.loc[data['description'].str.contains('programmed cell death 1 ', case=False)]

# Resume
print('Final programmed cell death 1 :')
print(len(data3.external_gene_name.value_counts()))

# Genes related with renal cell carcinoma
data4 = data.loc[data['description'].str.contains('renal cell carcinoma', case=False)]

# Resume
print('Final renal cell carcinoma:')
print(len(data4.external_gene_name.value_counts()))

# Genes related with prostate tumor
data5 = data.loc[data['description'].str.contains('prostate', case=False)]

# Resume
print('Final prostate:')
print(len(data5.external_gene_name.value_counts()))

# Genes related with bladder tumor
data6 = data.loc[data['description'].str.contains('bladder', case=False)]

# Resume
print('Final bladder:')
print(len(data6.external_gene_name.value_counts()))

# Genes related with wilms tumor
data7 = data.loc[data['description'].str.contains('wilms', case=False)]

# Resume
print('Final wilms:')
print(len(data7.external_gene_name.value_counts()))

data2 = pd.concat([data2,data3, data4, data5, data6, data7], axis=0)

data2.to_excel('GNN model\Data\Programmed cell death protein\Programmed cell death protein.xlsx', index=False)
