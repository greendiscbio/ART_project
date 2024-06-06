import pandas as pd
path ='Data/Clinical_data_of_patients_with_RNA_study.csv'
data = pd.read_csv(path, delimiter=';')
data = data[data['Arm'] == 'NIVOLUMAB']
data = data[['PFS', 'MSKCC', 'ImmunoPhenotype']]
data.MSKCC = data.MSKCC.astype("category").cat.codes
data.ImmunoPhenotype = data.ImmunoPhenotype.astype("category").cat.codes

print('Matriz de correlación entre PFS y MSKSCC:')
print(data.corr())


