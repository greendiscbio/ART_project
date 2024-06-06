
############################### METHODOLOGY REPRODUCTION OF #############################
################ "Convolutional neural network approach to lung cancer ##################
################ classification integrating protein interaction network #################
################               and gene expression profiles"            #################
#########################################################################################

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import networkx as nx
from sklearn.preprocessing import MinMaxScaler
import os
os.chdir("2023 Work/Graph networks/Biogrid")

# # print("################## Create eigenvalues and eigenvectors from edgelist ##################")

# # Get the graph of the biggest component
initial_G = nx.read_edgelist("../Data/HVGS/Spectral_clustering/HVGS/Biogrid/biogrid_found_genes_level_2_hvgs.edgelist")  
G = nx.Graph(initial_G)
print(G)
Gcc = sorted(nx.connected_components(G), key=len, reverse=True)
G = G.subgraph(Gcc[0])
print("Maximum component:")
print(G)
print()
# print("Adjacency Matrix:")
# A = nx.adjacency_matrix(G).todense()
# print(A)

# # diagonal matrix
# D = np.diag(A.sum(axis=1))

# # graph laplacian
# L = D-A

# # eigenvalues and eigenvectors
# vals, vecs = np.linalg.eig(L)

# # We only need the real part not the imaginary one (+0j)
# vals = vals.real
# vecs = vecs.real

# # sort these based on the eigenvalues
# vecs = vecs[:,np.argsort(vals)]
# vals = vals[np.argsort(vals)]

# # create data frame with the second and the third minimum eigenvalue
# df = pd.DataFrame(vecs.T, 
#                    columns=vals)
# df.to_csv("../Data/HVGS/Spectral_clustering/Data/HVGS/Spectral_clustering/Data/biggest_comp_eigenvalues_hvgs.csv")


print("################## Select the representative eigenvalues and eigenvectors ##################")

df = pd.read_csv("Data/biggest_comp_eigenvalues_hvgs.csv")
df = df.drop(columns=['Unnamed: 0'])
# print(df)

# Two first eigenvalues of the biggest component
df2 = df.iloc[:,3645:3647]
print('Eigenvalues: '+ str(df2.columns))
print(df2)

print("################## Create 2D representation ##################")

for i in range(len(df2)):
  plt.plot(df2.iloc[:,1][i], df2.iloc[:,0][i], marker="x",  markersize=7, color="green")
#plt.show()
plt.savefig("images/2d_representation_hvgs.png", format = 'png')
plt.cla()

print("################## Create 2D representation in the positive area ##################")

df2.iloc[:,1] = df2.iloc[:,1]+abs(df2.iloc[:,1].min())
df2.iloc[:,0] = df2.iloc[:,0]+abs(df2.iloc[:,0].min())
# print(df2)

for i in range(len(df2)):
  plt.plot(df2.iloc[:,1][i], df2.iloc[:,0][i], marker="x",  markersize=7, color="green")
#plt.show()
plt.savefig("../Data/HVGS/Spectral_clustering/images/2d_representation_positive_area_hvgs.png", format = 'png')
plt.cla()


# print("################## Reduce  space ##################")

## X axis
# df2 = df2[df2.iloc[:,1] < 0.35]
# df2 = df2[df2.iloc[:,1] > 0.15]
# ## Y axis
# df2 = df2[df2.iloc[:,0] < 0.20]
# df2 = df2[df2.iloc[:,0] > 0.0]
# # print(df2)

# # X axis
# df2 = df2[df2.iloc[:,1] < 0.3]
# df2 = df2[df2.iloc[:,1] > 0.225]
# ## Y axis
# df2 = df2[df2.iloc[:,0] < 0.11]
# df2 = df2[df2.iloc[:,0] > 0.08]

## More focus
df2 = df2[df2.iloc[:,1] < 0.275]
df2 = df2[df2.iloc[:,1] > 0.265]

## More focus
df2 = df2[df2.iloc[:,0] < 0.099]
df2 = df2[df2.iloc[:,0] > 0.095]

print(df2)
print(len(df2))
for i in df2.index:
  plt.plot(df2.iloc[:,1][i], df2.iloc[:,0][i], marker="x",  markersize=7, color="green")
#plt.show()

plt.savefig("../Data/HVGS/Spectral_clustering/images/2d_representation_positive_area_reduced_hvgs.png", format = 'png')
plt.cla()

df2 = df2*10000
df2 = df2.astype("int64")
print(df2)

print("#################### Create Heat map to see gene coord repetition #############################")

# Dimensionality reduction
df2.iloc[:,1] = df2.iloc[:,1]-abs(df2.iloc[:,1].min())
df2.iloc[:,0] = df2.iloc[:,0]-abs(df2.iloc[:,0].min())

heat = df2.groupby(df2.columns.tolist()).size().reset_index().rename(columns={0:'records'})

print('Maximum hit of the matrix: ' + str(heat["records"].max()))
print('Maximum value of the first eignevalue: ' + str(heat["0.16691354591782925"].max()))
print('Maximum value of the second eigenvalue: ' + str(heat["0.21990747247141856"].max()))
print('Minimum hit of the matrix: ' + str(heat["records"].min()))
print('Minimum value of the first eigenvalue: ' + str(heat["0.16691354591782925"].min()))
print('Minimum value of the second eigenvalue: ' + str(heat["0.21990747247141856"].min()))

print(df2)
print(heat)
df2.to_csv("../Data/HVGS/Spectral_clustering/Data/matrix_coords_zoom.csv")
heat.to_csv("../Data/HVGS/Spectral_clustering/Data/heat_coords_distribution_zoom.csv")

heats = np.ones((df2.iloc[:,1].max()+1, df2.iloc[:,0].max()+1))
heats[heat[df2.columns[1]], heat[df2.columns[0]]] = heat["records"]
print(heats.shape)

# import seaborn as sns
# sns.heatmap(heats, cmap ="coolwarm")
# #plt.show()
# plt.savefig("../Data/HVGS/Spectral_clustering/Data/heat_map.png", format = 'png')
# plt.cla()

path ='../../Data/Preprocessed_data/clinic_and_RNA_data_raw_NIVOLUMAB.csv'
data = pd.read_csv(path)
print(data)
X = []
Y = []

for patient in range(len(data)): # Recorro pacientes
  matrix = np.zeros((df2.iloc[:,1].max()+1, df2.iloc[:,0].max()+1))
  print(patient)
  for i in (df2.index):
    # print(df2.iloc[:,1][i])
    # print(df2.iloc[:,0][i])
    # print(list(G.nodes)[df2.index[i]])
    # print(data[list(G.nodes)[df2.index[i]]][patient])
    matrix[df2.iloc[:,1][i], df2.iloc[:,0][i]] = matrix[df2.iloc[:,1][i], df2.iloc[:,0][i]] + data[list(G.nodes)[i]][patient]
  final_matrix = matrix/heats
  scaler = MinMaxScaler()
  scaler.fit(final_matrix)
  final_matrix = scaler.transform(final_matrix)
  final_matrix=pd.DataFrame(final_matrix)
  final_matrix.to_csv('Data/Patients_zoom/final_matrix_'+str(patient)+'.csv')
  X.append(final_matrix)
  if data.PFS[patient] < 3:
    Y.append(0)
  else:
    Y.append(1)
print(X[0])