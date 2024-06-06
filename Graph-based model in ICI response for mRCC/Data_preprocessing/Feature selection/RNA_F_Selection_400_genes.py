import sklearn
#from sklearn.feature_selection import SequentialFeatureSelector as SFS
import numpy as np
#from sklearn.linear_model import LinearRegression
import pandas as pd
# Sequential Forward Selection(sfs)
#sfs = SFS(LinearRegression(),
 #         n_features_to_select=15,
#         direction = "forward",
#          scoring = 'f1')

from sklearn.neighbors import KNeighborsClassifier

from mlxtend.feature_selection import SequentialFeatureSelector as SFS
knn = KNeighborsClassifier(n_neighbors=4)

sfs = SFS(knn,
           k_features=400,
           forward=True,
           floating=False,
           verbose=1,
           scoring='accuracy',
           cv=0)
path ='RNA.csv'
data = pd.read_csv(path)
data = data.set_index('gene_name')
data = data.T
X = data.iloc[:,1:43894]

Y = []
for i in range (len(data)):
    if data.Target[i]=='R':
        Y.append(1)
    else:
        Y.append(0)


sfs.fit(X, Y)
print(sfs.k_feature_names_)
np.savetxt("400features.txt", sfs.k_feature_names_, fmt='%s')