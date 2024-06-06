import sklearn
import numpy as np
import pandas as pd
from sklearn.linear_model import LinearRegression
from mlxtend.feature_selection import SequentialFeatureSelector as SFS


sfs = SFS(LinearRegression(),
           k_features=300,
           forward=True,
           floating=True,
	       verbose = 1,
           cv=0)

# All patients
# path ='Data/clinic_and_RNA_data_raw.csv'
# NIVOLUMAB patients
path ='Data/clinic_and_RNA_data_raw_NIVOLUMAB.csv'

data = pd.read_csv(path)
print(data.head(5))

X = data.iloc[:,122:44017]
print(X)
Y = []
for i in range (len(data)):
    if data.PFS[i]<3:
        Y.append(0)
    else:
        Y.append(1)


sfs.fit(X, Y)
print(sfs.k_feature_names_)
# All patients
# np.savetxt("300_features.txt", sfs.k_feature_names_, fmt='%s')
# NIVOLUMAB patients
np.savetxt("300_features_NIVOLUMAB.txt", sfs.k_feature_names_, fmt='%s')

df_SFS_results = pd.DataFrame(sfs.subsets_).transpose()

# All patients
# df_SFS_results.to_csv('data300.csv')
# NIVOLUMAB patients
df_SFS_results.to_csv('data300_NIVOLUMAB.csv')