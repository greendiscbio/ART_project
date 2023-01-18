from sklearn.feature_selection import VarianceThreshold
import pandas as pd

path = "2023 Work/Data/Preprocessed_data/Gene Matrix/biogrid_included_genes.csv"
data_frame = pd.read_csv(path)
X = data_frame.iloc[:,1:17826 ]
print(X)

# Discarded Feature selection due to the low variance of the dataset.
# sel = VarianceThreshold(threshold=(.1 * (1 - .1)))
# print(sel.fit_transform(X).shape)

print("######################################")

from sklearn.feature_selection import SelectKBest
from sklearn.feature_selection import chi2

y = data_frame['Y']
print(X.shape)
num_features = input("How many features do you want to keep? ")
X_new = SelectKBest(chi2, k=int(num_features))
X_new.fit_transform(X, y)
cols = X_new.get_support(indices=True)
df = X.iloc[:,cols]
df['Y'] = y
df.to_csv('2023 Work/Graph networks/GNN+RN/Data/'+str(num_features)+'_most_relevant_genes_biogrid_found.csv')