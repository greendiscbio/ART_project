
############################### METHODOLOGY REPRODUCTION OF #############################
################ "Convolutional neural network approach to lung cancer ##################
################ classification integrating protein interaction network #################
################               and gene expression profiles"            #################
#########################################################################################

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import networkx as nx
import os
os.chdir("2023 Work/Graph networks/Biogrid")

# # print("################## Create eigenvalues and eigenvectors from edgelist ##################")

# # Get the graph of the biggest component
initial_G = nx.read_edgelist("Data/biogrid_found_genes_level_2_hvgs.edgelist")  
G = nx.Graph(initial_G)
print(G)
Gcc = sorted(nx.connected_components(G), key=len, reverse=True)
G = G.subgraph(Gcc[0])
print("Maximum component:")
print(G)
# # print()
# # print("Adjacency Matrix:")
# # A = nx.adjacency_matrix(G).todense()
# # print(A)

# # # diagonal matrix
# # D = np.diag(A.sum(axis=1))

# # # graph laplacian
# # L = D-A

# # # eigenvalues and eigenvectors
# # vals, vecs = np.linalg.eig(L)

# # # We only need the real part not the imaginary one (+0j)
# # vals = vals.real
# # vecs = vecs.real

# # # sort these based on the eigenvalues
# # vecs = vecs[:,np.argsort(vals)]
# # vals = vals[np.argsort(vals)]

# # # create data frame with the second and the third minimum eigenvalue
# # df = pd.DataFrame(vecs.T, 
# #                    columns=vals)
# # df.to_csv("Data/biggest_comp_eigenvalues_hvgs.csv")


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
plt.savefig("images/2d_representation_positive_area_hvgs.png", format = 'png')
plt.cla()


print("################## Reduce  space ##################")

## X axis
# df2 = df2[df2.iloc[:,1] < 0.35]
# df2 = df2[df2.iloc[:,1] > 0.15]
# ## Y axis
# df2 = df2[df2.iloc[:,0] < 0.20]
# df2 = df2[df2.iloc[:,0] > 0.0]
# # print(df2)

# X axis
df2 = df2[df2.iloc[:,1] < 0.3]
df2 = df2[df2.iloc[:,1] > 0.225]
## Y axis
df2 = df2[df2.iloc[:,0] < 0.11]
df2 = df2[df2.iloc[:,0] > 0.08]
print(df2)
print(len(df2))
for i in df2.index:
  plt.plot(df2.iloc[:,1][i], df2.iloc[:,0][i], marker="x",  markersize=7, color="green")
#plt.show()

plt.savefig("images/2d_representation_positive_area_reduced_hvgs.png", format = 'png')
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
# df2.to_csv("matrix_coords.csv")
# heat.to_csv("heat_coords_distribution.csv")

matrix = np.zeros((df2.iloc[:,1].max()+1, df2.iloc[:,0].max()+1))
matrix[heat[df2.columns[1]], heat[df2.columns[0]]] = heat["records"]
print(matrix)

import seaborn as sns
sns.heatmap(matrix, cmap ="coolwarm")
#plt.show()
plt.savefig("images/heat_map.png", format = 'png')
plt.cla()