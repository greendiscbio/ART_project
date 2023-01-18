import sklearn
import numpy as np
import pandas as pd
from sklearn.neighbors import KNeighborsClassifier
from mlxtend.feature_selection import SequentialFeatureSelector as SFS

knn = KNeighborsClassifier(n_neighbors=4)

sfs = SFS(knn,
           k_features=300,
           forward=True,
           floating=False,
           verbose=1,
           scoring='f1',
           cv=5)
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
np.savetxt("300_features_NIVOLUMAB.txt", sfs.k_feature_names_, fmt='%s')

