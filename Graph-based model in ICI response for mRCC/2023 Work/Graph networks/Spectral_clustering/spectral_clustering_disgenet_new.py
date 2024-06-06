
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
os.chdir("2023 Work/Graph networks/Spectral_clustering/")

# # print("################## Create eigenvalues and eigenvectors from edgelist ##################")
rna_total_file = pd.read_csv('../../Data/Preprocessed_data/clinic_and_RNA_data_raw_NIVOLUMAB.csv')
print(rna_total_file)
rna_genes = rna_total_file.columns.values[122:44015]
print('RNA_total genes: '+ str(len(rna_genes)))
# Get the graph of the biggest component
initial_G = nx.read_edgelist("Data/Patients Disgenet/biogrid_found_genes_level_3_disgenet.edgelist")  
G = nx.Graph(initial_G)
print(G)
not_included = []
for i in G.nodes():
  if i not in rna_genes:
    not_included.append(i)

print(not_included)
for i in not_included:
  G.remove_node(i)
print(G)
Gcc = sorted(nx.connected_components(G), key=len, reverse=True)
G = G.subgraph(Gcc[0])
print("Maximum component:")
print(G)

print("Adjacency Matrix:")
W = nx.adjacency_matrix(G)
print(W.todense())

# diagonal matrix
print('degree matrix:')
D = np.diag(np.sum(np.array(W.todense()), axis=1))
print(D)

# graph laplacian
print('laplacian matrix:')
L = D - W
print(L)

# eigenvalues and eigenvectors
e, v = np.linalg.eig(L)

# We only need the real part not the imaginary one (+0j)
e = e.real
v = v.real

# eigenvalues
print('eigenvalues:')
print(e)
# eigenvectors
print('eigenvectors:')
print(v)

# sort these based on the eigenvalues
v = v[:,np.argsort(e)]
e = e[np.argsort(e)]

# create data frame with the second and the third minimum eigenvalue
df = pd.DataFrame(v.T, columns=e)
print(df)
df.to_csv("Data/Patients Disgenet/biggest_comp_eigenvalues_3_disgenet.csv")


print("################## Select the representative eigenvalues and eigenvectors ##################")

df = pd.read_csv("Data/Patients Disgenet/biggest_comp_eigenvalues_3_disgenet.csv")
df = df.drop(columns=['Unnamed: 0'])

# Two first eigenvalues of the biggest component
print('Eigenvalues: '+ str(df.columns))

print("################## Create 2D representation ##################")

for i in range(len(df)):
  plt.plot(df.iloc[:,1][i], df.iloc[:,2][i], marker="x",  markersize=7, color="green")
#plt.show()
plt.savefig("Data/Patients Disgenet/images/2d_representation_disgenet_3.png", format = 'png')
plt.cla()

print("################## Create 2D representation in the positive area ##################")

df.iloc[:,1] = df.iloc[:,1]+abs(df.iloc[:,1].min())
df.iloc[:,2] = df.iloc[:,2]+abs(df.iloc[:,2].min())
print(df)

for i in range(len(df)):
  plt.plot(df.iloc[:,1][i], df.iloc[:,2][i], marker="x",  markersize=7, color="green")
# plt.show()
plt.plot(0.141,0.676, marker=".",  markersize=7, color = 'red')
plt.plot(0.148,0.672, marker=".",  markersize=7, color = 'red')
plt.plot(0.151,0.678, marker=".",  markersize=7, color = 'red')
plt.plot(0.204,0.687, marker=".",  markersize=7, color = 'red')
plt.plot(0.144,0.677, marker=".",  markersize=7, color = 'red')
plt.plot(0.145,0.676, marker=".",  markersize=7, color = 'red')
plt.plot(0.144,0.678, marker=".",  markersize=7, color = 'red')
plt.plot(0.142,0.678, marker=".",  markersize=7, color = 'red')
plt.plot(0.149,0.695, marker=".",  markersize=7, color = 'red')
plt.plot(0.145,0.679, marker=".",  markersize=7, color = 'red')

plt.title('Spectral clusteing result with most representative genes.')
plt.savefig("Data/Patients Disgenet/images/2d_representation_positive_area_disgenet_3.png", format = 'png')
plt.cla()


print("################## Reduce  space ##################")

# # X axis
df = df[df.iloc[:,1] < 0.205]
df = df[df.iloc[:,1] > 0.12]
## Y axis
df = df[df.iloc[:,2] < 0.7]
df = df[df.iloc[:,2] > 0.65]

print(df)
print(len(df))
for i in df.index:
  plt.plot(df.iloc[:,1][i], df.iloc[:,2][i], marker="x",  markersize=7, color="green")
#plt.show()
plt.plot(0.141,0.676, marker=".",  markersize=7, color = 'red')
plt.plot(0.148,0.672, marker=".",  markersize=7, color = 'red')
plt.plot(0.151,0.678, marker=".",  markersize=7, color = 'red')
plt.plot(0.204,0.687, marker=".",  markersize=7, color = 'red')
plt.plot(0.144,0.677, marker=".",  markersize=7, color = 'red')
plt.plot(0.145,0.676, marker=".",  markersize=7, color = 'red')
plt.plot(0.144,0.678, marker=".",  markersize=7, color = 'red')
plt.plot(0.142,0.678, marker=".",  markersize=7, color = 'red')
plt.plot(0.149,0.695, marker=".",  markersize=7, color = 'red')
plt.plot(0.145,0.679, marker=".",  markersize=7, color = 'red')

plt.title('Spectral clusteing result with most representative genes.')
plt.savefig("Data/Patients Disgenet/images/2d_representation_positive_area_reduced_disgenet_3.png", format = 'png')
# plt.cla()

df = df*1000
df = df.astype("int64")
print(df)

print("#################### Create Heat map to see gene coord repetition #############################")

# Dimensionality reduction
df.iloc[:,1] = df.iloc[:,1]-abs(df.iloc[:,1].min())
df.iloc[:,2] = df.iloc[:,2]-abs(df.iloc[:,2].min())
df = df.iloc[:,1:3]

heat = df.groupby(df.columns.tolist()).size().reset_index().rename(columns={0:'records'})

print('Maximum hit of the matrix: ' + str(heat["records"].max()))
print('Maximum value of the first eignevalue: ' + str(heat["0.0766620426658494"].max()))
print('Maximum value of the second eigenvalue: ' + str(heat["0.08043809987034627"].max()))
print('Minimum hit of the matrix: ' + str(heat["records"].min()))
print('Minimum value of the first eigenvalue: ' + str(heat["0.0766620426658494"].min()))
print('Minimum value of the second eigenvalue: ' + str(heat["0.08043809987034627"].min()))

print(df)
print(heat)
df.to_csv("Data/Patients Disgenet/matrix_coords_zoom.csv")
heat.to_csv("Data/Patients Disgenet/heat_coords_distribution_zoom.csv")

heats = np.ones((df.iloc[:,0].max()+1, df.iloc[:,1].max()+1))
heats[heat[df.columns[0]], heat[df.columns[1]]] = heat["records"]
print(heats.shape)

# # import seaborn as sns
# # sns.heatmap(heats, cmap ="coolwarm")
# # #plt.show()
# # plt.savefig("../Data/HVGS/Spectral_clustering/Data/heat_map.png", format = 'png')
# # plt.cla()

path ='../../Data/Preprocessed_data/clinic_and_RNA_data_raw_NIVOLUMAB.csv'
data = pd.read_csv(path)
print(data)
X = []
Y = []

for patient in range(len(data)): # Recorro pacientes
  matrix = np.zeros((df.iloc[:,0].max()+1, df.iloc[:,1].max()+1))
  print(patient)
  for i in (df.index):
    # print(df.iloc[:,1][i])
    # print(df.iloc[:,0][i])
    # print(i)
    # print(data[list(G.nodes)[df.index[i]]][patient])
    matrix[df.iloc[:,0][i], df.iloc[:,1][i]] = matrix[df.iloc[:,0][i], df.iloc[:,1][i]] + data[list(G.nodes)[i]][patient]
    with open('aa3.txt', 'a') as file:
      file.write(str(df.iloc[:,0][i])+', '+ str(df.iloc[:,1][i])+', '+ str(list(G.nodes)[i])+'\n')
  final_matrix = matrix/heats
  scaler = MinMaxScaler()
  scaler.fit(final_matrix)
  final_matrix = scaler.transform(final_matrix)
  final_matrix=pd.DataFrame(final_matrix)
  final_matrix.to_csv('Data/Patients Disgenet/patients/final_matrix_'+str(patient)+'.csv')
  X.append(final_matrix)
  if data.PFS[patient] < 3:
    Y.append(0)
  else:
    Y.append(1)
print(X[0])
# print(G.nodes)