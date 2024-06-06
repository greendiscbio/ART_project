import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import statsmodels.api as sm

path ='Data\Clinical_data_of_patients_with_RNA_study.csv'
data = pd.read_csv(path,  delimiter=';')

CM_009 = data[data.Cohort=="CM-009"]
CM_010 = data[data.Cohort=="CM-010"]
CM_025 = data[data.Cohort=="CM-025"]



result = [len(CM_009[CM_009.Arm=="NIVOLUMAB"]), len(CM_009[CM_009.Arm=="EVEROLIMUS"])]
result2 = [len(CM_010[CM_010.Arm=="NIVOLUMAB"]), len(CM_010[CM_010.Arm=="EVEROLIMUS"])]
result3 = [len(CM_025[CM_025.Arm=="NIVOLUMAB"]), len(CM_025[CM_025.Arm=="EVEROLIMUS"])]

names = ["Nivolumab", "Everolimus"]

fig, ((ax1, ax2, ax3)) = plt.subplots(1, 3)
fig.suptitle('Cohorts description')
ax1.pie(result, labels=names, autopct="%i %%")
ax2.pie(result2, labels=names, autopct="%i %%")
ax3.pie(result3, labels=names, autopct="%i %%")
ax1.set_title('Patients in cohort CM-009: '+ str(len(CM_009)))
ax2.set_title('Patients in cohort CM-010: '+ str(len(CM_010)))
ax3.set_title('Patients in cohort CM-025: '+ str(len(CM_025)))

plt.show()
