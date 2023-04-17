import pandas as pd
import numpy as np


path ='2023 Work/Data/Preprocessed_data/clinic_and_RNA_data_raw_NIVOLUMAB.csv'
data = pd.read_csv(path)
print('Init')
median_r_matrix = np.zeros((100,41))
median_nr_matrix = np.zeros((100,41))
cont_r = 0
cont_nr = 0
for i in range(180):
    patient = pd.read_csv('2023 Work/Graph networks/Spectral_clustering/Data/Patients_zoom/final_matrix_'+str(i)+'.csv')
    if data.PFS[i] < 3:
        median_nr_matrix = median_nr_matrix + patient.to_numpy()
        cont_nr = cont_nr + 1
    else:
        median_r_matrix = median_r_matrix + patient.to_numpy()
        cont_r = cont_r + 1
print(cont_r)
print(cont_nr)
print(median_r_matrix)
final_matrix=pd.DataFrame(median_r_matrix/cont_r)
print('Responder patients data:')
print('Max: '+str((median_r_matrix/cont_r).max()))
print('Min: '+str((median_r_matrix/cont_r).min()))

print('Non-Responder patients data:')
print('Max: '+str((median_nr_matrix/cont_nr).max()))
print('Min: '+str((median_nr_matrix/cont_nr).min()))
final_matrix.to_csv('2023 Work/Graph networks/Spectral_clustering/Data/Patients_zoom/R_patients.csv')
final_matrix=pd.DataFrame(median_nr_matrix/cont_nr)
final_matrix.to_csv('2023 Work/Graph networks/Spectral_clustering/Data/Patients_zoom/NR_patients.csv')

nr_patient = pd.read_csv('2023 Work/Graph networks/Spectral_clustering/Data/Patients_zoom/NR_patients.csv').to_numpy()
r_patient = pd.read_csv('2023 Work/Graph networks/Spectral_clustering/Data/Patients_zoom/R_patients.csv').to_numpy()

deviation = nr_patient-r_patient
print('Deviation: ')
print('Max: '+str(deviation.max()))
print('Min: '+str(deviation.min()))
final_matrix=pd.DataFrame(nr_patient-r_patient)

final_matrix.to_csv('2023 Work/Graph networks/Spectral_clustering/Data/Patients_zoom/std_patients.csv')
