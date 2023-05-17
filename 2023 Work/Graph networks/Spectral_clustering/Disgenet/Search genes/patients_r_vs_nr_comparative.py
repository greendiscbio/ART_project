import pandas as pd
import numpy as np

path ='2023 Work/Data/Preprocessed_data/clinic_and_RNA_data_raw_NIVOLUMAB.csv'
data = pd.read_csv(path)

median_r_matrix = np.zeros((1062,353))
median_nr_matrix = np.zeros((1062,353))

cont_r = 0
cont_nr = 0

for i in range(181):
    print('Patient '+str(i))
    patient = pd.read_csv('2023 Work/Graph networks/Spectral_clustering/Data/Patients Disgenet/patients/final_matrix_'+str(i)+'.csv')
    patient = patient.drop(['Unnamed: 0'], axis=1)
    # print(patient)
    if data.PFS[i] < 3:
        median_nr_matrix = median_nr_matrix + patient.to_numpy()
        cont_nr = cont_nr + 1
    else:
        median_r_matrix = median_r_matrix + patient.to_numpy()
        cont_r = cont_r + 1

final_matrix=pd.DataFrame(median_r_matrix/cont_r)
# .iloc[:,1:435]
print('Responder patients data: ' + str(cont_r)+ ' patients')
print('Max: '+str((final_matrix.to_numpy()).max()))
print('Min: '+str((final_matrix.to_numpy()).min()))
final_matrix.to_csv('2023 Work/Graph networks/Spectral_clustering/Data/Patients Disgenet/patients/R_patients.csv')

final_matrix=pd.DataFrame(median_nr_matrix/cont_nr)
# .iloc[:,1:435]
print('Non-Responder patients data: '+ str(cont_nr)+ ' patients')
print('Max: '+str((final_matrix.to_numpy()).max()))
print('Min: '+str((final_matrix.to_numpy()).min()))
final_matrix.to_csv('2023 Work/Graph networks/Spectral_clustering/Data/Patients Disgenet/patients/NR_patients.csv')

nr_patient = pd.read_csv('2023 Work/Graph networks/Spectral_clustering/Data/Patients Disgenet/patients/NR_patients.csv').to_numpy()
r_patient = pd.read_csv('2023 Work/Graph networks/Spectral_clustering/Data/Patients Disgenet/patients/R_patients.csv').to_numpy()

deviation = nr_patient-r_patient
print('Deviation: ')
print('Max deviation: '+str(deviation.max()) + ' found at '+str(np.unravel_index(np.argmax(deviation, axis=None), deviation.shape)))
print('Min: '+str(deviation.min())+ ' found at '+str(np.unravel_index(np.argmin(deviation, axis=None), deviation.shape)))
final_matrix=pd.DataFrame(nr_patient-r_patient)

final_matrix.to_excel('2023 Work/Graph networks/Spectral_clustering/Data/Patients Disgenet/patients/std_patients.xlsx')