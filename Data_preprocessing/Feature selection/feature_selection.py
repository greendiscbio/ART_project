import numpy as np
from sklearn.linear_model import LinearRegression
import pandas as pd
from mlxtend.feature_selection import SequentialFeatureSelector as SFS
from sklearn.model_selection import StratifiedKFold # import KFold

kf=StratifiedKFold(n_splits=5, random_state=None, shuffle=False)

path ='../clinic_and_RNA_data_raw_NIVOLUMAB.csv'
data = pd.read_csv(path)
print(data.head())

data['PFS_label'] = pd.cut(data['PFS'], bins=[0, 3, np.inf], labels=['NR', 'R'])
Y = data['PFS_label'].replace({'NR':0, 'R':1}).values
X = data.iloc[:,122:44015]

sfs = SFS(LinearRegression(),
            k_features=150,
            forward=True,
            floating=True,
            verbose = 1,
            cv=5)
for train_index, val_index in kf.split(X, Y):
    i = i + 1
    print("TRAIN: ", train_index, "VAL:", val_index)

    X_train, X_val = X.iloc[train_index], X.iloc[val_index]
    Y_train, Y_val =  np.array(Y)[train_index],  np.array(Y)[val_index]
    
    sfs.fit(X_train, Y_train)
    print(sfs.k_feature_names_)
    df_SFS_results = pd.DataFrame(sfs.subsets_).transpose()
    df_SFS_results.to_csv('150_features_Selection.csv')
