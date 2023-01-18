############################### METHODOLOGY REPRODUCTION OF #############################
# ############## "Convolutional neural network approach to lung cancer ##################
# ############## classification integrating protein interaction network# ################
################               and gene expression profiles"             ################
# #######################################################################################

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import networkx as nx

initial_G = nx.read_edgelist("2023 Work/Graph networks/Biogrid/Data/biogrid_minimum.edgelist")  
G = nx.Graph(initial_G)
A = nx.adjacency_matrix(G).todense()

print("Adjacency Matrix:")
print(A)

# diagonal matrix
D = np.diag(A.sum(axis=1))

# graph laplacian
L = D-A

# eigenvalues and eigenvectors
vals, vecs = np.linalg.eig(L)

# create data frame with the second and the third minimum eigenvalue
df = pd.DataFrame(vecs.T, 
                   columns=vals)
df.to_csv("2023 Work\Graph networks\Biogrid\Data\eigenvalues_RNA_30.csv")
print(df)
print(vals)
first_eigen_val = df.columns[np.where(df.columns == np.amin(df.columns))[0]][0]
df = df.drop([first_eigen_val], axis=1)

second_eigen_val = df.columns[np.where(df.columns == np.amin(df.columns))[0]][0]
df2 = pd.DataFrame(df[second_eigen_val], columns=[second_eigen_val])
df = df.drop([second_eigen_val], axis=1)

third_eigen_val = df.columns[np.where(df.columns == np.amin(df.columns))[0]][0]
col = df[third_eigen_val]

df2[third_eigen_val] =  col

print("################## Create 2D representation ###############################")
print(df2)
print(np.abs(np.round(df2[second_eigen_val])).max()+1)
print(np.abs(np.round(df2[third_eigen_val])).max()+1)

for i in range(len(df2)):
  plt.plot(df2[third_eigen_val][i], df2[second_eigen_val][i], marker="x",  markersize=7, color="green")

plt.show()

print("#################### Create 2D representation in the positive area #############################")

print(df2[third_eigen_val].min())
df2[third_eigen_val] = df2[third_eigen_val]+abs(df2[third_eigen_val].min())
df2[second_eigen_val] = df2[second_eigen_val]+abs(df2[second_eigen_val].min())
print(df2)
for i in range(len(df2)):
  plt.plot(df2[third_eigen_val][i], df2[second_eigen_val][i], marker="x",  markersize=7, color="green")
plt.show()

df2 = df2*100
df2 = df2.astype("int32")
print(df2)
matrix = np.zeros((df2[third_eigen_val].max()+1, df2[second_eigen_val].max()+1))
print(matrix.size)
# df1 = df2[df2.duplicated(keep=False)]

# print (df1)

print("#################### Create Heat map to see gene coord repetition #############################")

heat = df2.groupby(df2.columns.tolist()).size().reset_index().rename(columns={0:'records'})
print(heat)

matrix[heat[third_eigen_val], heat[second_eigen_val]] = heat["records"]
print(matrix)
import seaborn as sns
sns.heatmap(matrix)
plt.show()

######################## GENE MAPPING ########################


matrix = np.zeros((100, 100))
path ='2023 Work\Data\Preprocessed_data\clinic_and_RNA_data_raw.csv'
data = pd.read_csv(path)
# print(data)
X = []
Y = []
for patient in range(len(data)): # Recorro pacientes
  print(patient)
  for i in range(len(df2)):
    # print(df2.iloc[i])
    # print(df2[third_eigen_val][i])
    # print(df2[second_eigen_val][i])
    # print(list(G.nodes)[i])
    # print(data[list(G.nodes)[i]][patient])
    matrix[df2[third_eigen_val][i], df2[second_eigen_val][i]] = data[list(G.nodes)[i]][patient]
  X.append(matrix)
  if data.PFS[i] < 3:
    Y.append(0)
  else:
    Y.append(1)
print(X[0])

#################################### CNN #####################################

# from sklearn.model_selection import train_test_split 
# train_data, test_data, Y_train, Y_test = train_test_split(X, Y, test_size=0.1, stratify=Y)

# from sklearn.model_selection import StratifiedKFold # import KFold
# kf=StratifiedKFold(n_splits=10, random_state=None, shuffle=False)

# val_avg = []
# test_avg = []
# test_f1_score = []
# for train_index, val_index in kf.split(train_data, Y_train):
#     train_dataset=[]
#     val_dataset=[]
#     print("TRAIN: ", train_index, "TEST:", val_index)
#     for i in train_index:
#         train_dataset.append(train_data[i])
#     for i in val_index:
#         val_dataset.append(train_data[i])

#     print(len(train_dataset))
#     print(len(val_dataset))

# import tensorflow as tf

# from tensorflow.keras import datasets, layers, models

# model = models.Sequential()
# model.add(layers.Conv2D(512, (5, 5), activation='relu', input_shape=(100, 100, 1)))
# model.add(layers.MaxPooling2D((2, 2)))
# model.add(layers.Conv2D(256, (3, 3), activation='relu'))
# model.add(layers.MaxPooling2D((2, 2)))
# model.add(layers.Conv2D(128, (3, 3), activation='relu'))
# model.add(layers.MaxPooling2D((2, 2)))
# model.add(layers.Dense(256))
# model.add(layers.Dense(256))
# model.add(layers.Dense(128))

# print(model.output_shape)