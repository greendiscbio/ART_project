from sklearn.cluster import KMeans
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
# our adjacency matrix


# Adjacency Matrix:
# A = np.array([
# [0, 1, 1, 0, 0, 1, 0, 0, 1, 1,],
#  [1, 0, 1, 0, 0, 0, 0, 0, 0, 0],
#  [1, 1, 0, 0, 0, 0, 0, 0, 0, 0],
#  [0, 0, 0, 0, 1, 1, 0, 0, 0, 0],
#  [0, 0, 0, 1, 0, 1, 0, 0, 0, 0],
#  [1, 0, 0, 1, 1, 0, 1, 1, 0, 0],
#  [0, 0, 0, 0, 0, 1, 0, 1, 0, 0],
#  [0, 0, 0, 0, 0, 1, 1, 0, 0, 0],
#  [1, 0, 0, 0, 0, 0, 0, 0, 0, 1],
#  [1, 0, 0, 0, 0, 0, 0, 0, 1, 0]])
A = np.array([
  [0, 1, 1, 0, 0, 0, 0, 0, 1, 1],
  [1, 0, 1, 0, 0, 0, 0, 0, 0, 0],
  [1, 1, 0, 0, 0, 0, 0, 0, 0, 0],
  [0, 0, 0, 0, 1, 1, 0, 0, 0, 0],
  [0, 0, 0, 1, 0, 1, 0, 0, 0, 0],
  [0, 0, 0, 1, 1, 0, 1, 1, 0, 0],
  [0, 0, 0, 0, 0, 1, 0, 1, 0, 0],
  [0, 0, 0, 0, 0, 1, 1, 0, 0, 0],
  [1, 0, 0, 0, 0, 0, 0, 0, 0, 1],
  [1, 0, 0, 0, 0, 0, 0, 0, 1, 0]])

print("Adjacency Matrix:")
print(A)

# diagonal matrix
D = np.diag(A.sum(axis=1))

# graph laplacian
L = D-A

# eigenvalues and eigenvectors
vals, vecs = np.linalg.eig(L)

# sort these based on the eigenvalues
vecs = vecs[:,np.argsort(vals)]
vals = vals[np.argsort(vals)]

print(vals)
# print(vecs)
df = pd.DataFrame(vecs.T, 
                   columns=vals)
print(df)
df2 = df.iloc[:,2:4]
print(df2)
print(df2.columns)
# print(df2.loc(0))
# plt.scatter(df2[df2.columns[0]],df2[df2.columns[1]])
# for i in range(len(df2)):
#     x = df2[df2.columns[0]][i]
#     y = df2[df2.columns[1]][i]
#     plt.annotate(i, (x,y), textcoords="offset points", xytext=(0,10))
#     print(i)
#     print(x,y)
#   plt.plot(df2[df2.columns[1]][i], df2[df2.columns[0]][i], marker="x",  markersize=7, color="green")
#   print(df2[df2.columns[0]][i], df2[df2.columns[1]][i])
# plt.show()

# kmeans on first three vectors with nonzero eigenvalues
kmeans = KMeans(n_clusters=2)
kmeans.fit(vecs[:,2:4])
colors = kmeans.labels_

print("Clusters:", colors)