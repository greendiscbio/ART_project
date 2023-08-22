################### Search for included genes in Disgenet, RNA_total and Biogrid ###################

import pandas as pd
import numpy as np
disgenet_genes = []
disgenet_file = pd.read_excel('2023 Work/Data/Disgenet_diasease_genes.xlsx')
print(disgenet_file)
disgenet_file=disgenet_file[disgenet_file['Disease'].isin( ['Conventional (Clear Cell) Renal Cell Carcinoma','Clear-cell metastatic renal cell carcinoma', 'Hereditary clear cell renal cell carcinoma', 'Non-Hereditary Clear Cell Renal Cell Carcinoma'])]
disgenet_genes = disgenet_file['Gene'].unique()
print('Disgenet genes: '+str(len(disgenet_genes)))

rna_total_genes = []
rna_total_file = pd.read_csv('2023 Work/Data/Preprocessed_data/clinic_and_RNA_data_raw_NIVOLUMAB.csv')
print(rna_total_file)
rna_total_genes = rna_total_file.columns.values[122:44015]
print('RNA_total genes: '+ str(len(rna_total_genes)))

comparation = []


for item in disgenet_genes:
  if item in rna_total_genes: #and item in biogrid_genes:
    comparation.append(item)
    
print(len(comparation))
print('Disgenet included genes: '+ str(len(comparation)))
exit()
with open("2023 Work/Graph networks/Graph2Vec/ppi_creation/Data/disgenet_found_genes_in_rna_total_Nivolumab.txt", 'w') as file:
    for c in comparation:
        file.write(c+'\n')


rna_total_genes = []
rna_total_file = pd.read_csv('2023 Work/Data/Preprocessed_data/clinic_and_RNA_data_raw_Ave_Axi.csv')
print(rna_total_file)
rna_total_genes = rna_total_file.columns.values[122:44015]
print('RNA_total genes: '+ str(len(rna_total_genes)))

comparation = []


for item in disgenet_genes:
  if item in rna_total_genes: #and item in biogrid_genes:
    comparation.append(item)
    
print(len(comparation))
print('Disgenet included genes: '+ str(len(comparation)))
with open("2023 Work/Graph networks/Graph2Vec/ppi_creation/Data/disgenet_found_genes_in_rna_total_Ave_Axi.txt", 'w') as file:
    for c in comparation:
        file.write(c+'\n')


################### Search for deviation in Disgenet genes ###################

# rna_total_file = rna_total_file.loc[:, comparation]
# rna_total_file.to_csv('RNA_Disgenet_matrix.csv')



# median_r_matrix = np.zeros((2144, 1))
# median_nr_matrix = np.zeros((2144, 1))
# cont_r = 0
# cont_nr = 0

# for patient in rna_total_file.iloc:
#    if patient.PFS < 3:
#         median_nr_matrix = median_nr_matrix + patient.to_numpy()
#         cont_nr = cont_nr + 1
#    else:
#         median_r_matrix = median_r_matrix + patient.to_numpy()
#         cont_r = cont_r + 1
# print(median_r_matrix)
# final_matrix=pd.DataFrame(median_r_matrix/cont_r).iloc[0].T
# print(final_matrix)
# print('Responder patients data: ' + str(cont_r)+ ' patients')
# print('Max: '+str((final_matrix.to_numpy()).max()))
# print('Min: '+str((final_matrix.to_numpy()).min()))
# final_matrix.to_csv('R_patients.csv')

# final_matrix=pd.DataFrame(median_nr_matrix/cont_nr).iloc[0].T
# print('Non-Responder patients data: '+ str(cont_nr)+ ' patients')
# print('Max: '+str((final_matrix.to_numpy()).max()))
# print('Min: '+str((final_matrix.to_numpy()).min()))
# final_matrix.to_csv('NR_patients.csv')

# nr_patient = pd.read_csv('NR_patients.csv').to_numpy()
# r_patient = pd.read_csv('R_patients.csv').to_numpy()

# deviation = nr_patient-r_patient
# print('Deviation: ')
# print('Max deviation: '+str(deviation.max()) + ' found at '+str(np.unravel_index(np.argmax(deviation, axis=None), deviation.shape)))
# print('Min: '+str(deviation.min())+ ' found at '+str(np.unravel_index(np.argmin(deviation, axis=None), deviation.shape)))
# final_matrix=pd.DataFrame(nr_patient-r_patient)

# final_matrix = final_matrix.T
# final_matrix.columns = comparation
# final_matrix = final_matrix.T
# final_matrix.to_excel('std_patients.xlsx')