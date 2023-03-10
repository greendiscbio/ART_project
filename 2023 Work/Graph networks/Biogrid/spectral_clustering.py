
############################### METHODOLOGY REPRODUCTION OF #############################
################ "Convolutional neural network approach to lung cancer ##################
################ classification integrating protein interaction network #################
################               and gene expression profiles"             ################
#########################################################################################

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import networkx as nx
import os
os.chdir("2023 Work/Graph networks/Biogrid")
# print("################## Create eigenvalues and eigenvectors from edgelist ##################")

# # Get the graph of the biggest component
# initial_G = nx.read_edgelist("Data/biogrid_found_genes_level_2_anova.edgelist")  
# G = nx.Graph(initial_G)
# A = nx.adjacency_matrix(G).todense()
# Gcc = sorted(nx.connected_components(G), key=len, reverse=True)
# G = G.subgraph(Gcc[0])
# print("Maximum component:")
# print(G)
# print()
# print("Adjacency Matrix:")
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
# df.to_csv("Data/biggest_comp_eigenvalues_anova.csv")


print("################## Select the representative eigenvalues and eigenvectors ##################")

df = pd.read_csv("Data/biggest_comp_eigenvalues_anova.csv")
df = df.drop(columns=['Unnamed: 0'])
# print(df)
# print(df.columns)
#print(vals)

#ANOVA
df2 = df.iloc[:,3769:3771]
#ANOVA Biggest comp
# df2 = df.iloc[:,3770:3772]

# # print(df2.columns)

# second_eigenvalue = df2.iloc[:,0]
# third_eigenvalue = df2.iloc[:,1]


print("################## Create 2D representation ##################")

# print(np.abs(np.round(df2.iloc[:,0])).max()+1)
# print(np.abs(np.round(df2.iloc[:,1])).max()+1)

for i in range(len(df2)):
  plt.plot(df2.iloc[:,1][i], df2.iloc[:,0][i], marker="x",  markersize=7, color="green")
#plt.show()
plt.savefig("images/2d_representation.png", format = 'png')
plt.cla()


print("################## Create 2D representation in the positive area ##################")

# print(df2.iloc[:,1].min())
df2.iloc[:,1] = df2.iloc[:,1]+abs(df2.iloc[:,1].min())
df2.iloc[:,0] = df2.iloc[:,0]+abs(df2.iloc[:,0].min())
# print(df2)
for i in range(len(df2)):
  plt.plot(df2.iloc[:,1][i], df2.iloc[:,0][i], marker="x",  markersize=7, color="green")
#plt.show()
plt.savefig("images/2d_representation_positive_area.png", format = 'png')
plt.cla()


print("################## Reduce  space ##################")

# print(df2)
## X axis
df2 = df2[df2.iloc[:,1] < 0.4]
# print(len(df2))
df2 = df2[df2.iloc[:,1] > 0.312]
# print(len(df2))
## Y axis
df2 = df2[df2.iloc[:,0] < 0.6]
# print(len(df2))
df2 = df2[df2.iloc[:,0] > 0.4]
# print(len(df2))
# print(df2)
df2 = df2.reset_index()
# df2 = df2.set_index("index")
df2.to_csv("Data/aa.csv")
# print(df2.iloc[:,1])
# print(df2.iloc[:,0])
for i in range(len(df2)):
  # print(df2.iloc[:,2][i], df2.iloc[:,1][i])
  plt.plot(df2.iloc[:,2][i], df2.iloc[:,1][i], marker="x",  markersize=7, color="green")
#plt.show()
plt.savefig("images/2d_representation_positive_area_reduced.png", format = 'png')
plt.cla()
print("aaaaaaaaaaaaaaaaaaaa")
print(df2)
df2 = df2*10000000000000000000
df2 = df2.astype("int64")
print(df2)
df2.to_csv("df2.csv")
# matrix = np.zeros((df2.iloc[:,2].max()+1, df2.iloc[:,1].max()+1))
# print(matrix.size)
# df1 = df2[df2.duplicated(keep=False)]
# print (df1)


print("#################### Create Heat map to see gene coord repetition #############################")

df2 = df2.drop(columns=['index'])
heat = df2.groupby(df2.columns.tolist()).size().reset_index().rename(columns={0:'records'})
print(heat)
heat.to_csv("heats.csv")
df2.to_csv("df2.csv")

# # TRY
heat = pd.read_csv("heats.csv")
df2 = pd.read_csv("df2.csv")

print(df2.iloc[:,1].min())
df2.iloc[:,2] = df2.iloc[:,2]-abs(df2.iloc[:,2].min())
df2.iloc[:,1] = df2.iloc[:,1]-abs(df2.iloc[:,1].min())
print(df2)
df2 = df2.drop(columns=['Unnamed: 0'])
heat = df2.groupby(df2.columns.tolist()).size().reset_index().rename(columns={0:'records'})

# # for i in range(len(df2)):
# #   plt.plot(df2.iloc[:,2][i], df2.iloc[:,1][i], marker="x",  markersize=7, color="green")
# # #plt.show()
# # plt.savefig("images/2d_representation_positive_area_after_all.png", format = 'png')
# # plt.cla()
# df2 = df2.drop(columns=['Unnamed: 0'])
# heat = df2.groupby(df2.columns.tolist()).size().reset_index().rename(columns={0:'records'})
# print(heat["records"])
print(heat["records"].max())
print(heat["0.23856184961504265"].max())
print(heat["0.2449362200197798"].max())
print(heat["records"].min())
print(heat["0.23856184961504265"].min())
print(heat["0.2449362200197798"].min())

heat.to_csv("heats.csv")
df2.to_csv("df2.csv")
heat = pd.read_csv("heats.csv")
df2 = pd.read_csv("df2.csv")
# matrix = np.zeros((df2.iloc[:,2].max()+1, df2.iloc[:,1].max()+1))
# matrix[heat[df2.columns[2]], heat[df2.columns[1]]] = heat["records"]
# print(matrix)

# import seaborn as sns
# sns.heatmap(matrix, cmap ="coolwarm")
# #plt.show()
# plt.savefig("images/heat_map.png", format = 'png')
# plt.cla()




######################## GENE MAPPING ########################

# print("####################### GENE MAPPING ###################")
# matrix = np.zeros((10000000000000000000, 10000000000000000000))
# path ='../Data/clinic_and_RNA_data_raw.csv'
# data = pd.read_csv(path)
# # print(data)
# X = []
# Y = []
# for patient in range(len(data)): # Recorro pacientes
#   print(patient)
#   for i in range(len(df2)):
#     # print(df2.iloc[i])
#     # print(df2[third_eigen_val][i])
#     # print(df2[second_eigen_val][i])
#     # print(list(G.nodes)[i])
#     # print(data[list(G.nodes)[i]][patient])
#     matrix[df2.iloc[:,1], df2.iloc[:,0]] = data[list(G.nodes)[i]][patient]
#   X.append(matrix)
#   if data.PFS[patient] < 3:
#     Y.append(0)
#   else:
#     Y.append(1)
# print(X[0])


# vals = []
# for i in range(0, 1019938561312778240):
#     vals.append(i)

# df = pd.DataFrame(columns = vals)
# df['x'] = vals
# df = df.set_index('x')
# print(df)
# matrix = pd.DataFrame(columns=columns)
# matrix.set_index()
# print(matrix)
# matrix.to_csv("matrix.csv")