import numpy as np
import matplotlib.pyplot as plt
import statistics
data = []
f = open("RNA_pfs.txt", "r")
for linea in f:
   l = float(linea[:-2])
   data.append(l)
f.close()

data2 = []
f = open("WES_pfs.txt", "r")
for linea in f:
   l = float(linea[:-2])
   data2.append(l)
f.close()

data3 = []
f = open("RNA_WES_pfs.txt", "r")
for linea in f:
   l = float(linea[:-2])
   data3.append(l)
f.close()

data4 = []
f = open("Clinico_pfs.txt", "r")
for linea in f:
   l = float(linea[:-2])
   data4.append(l)
f.close()

inter = np.arange(70)

fig, ((ax1, ax4), (ax3, ax2)) = plt.subplots(2, 2)
fig.suptitle('PFS en distintas tablas')
ax1.hist(x = data, bins = inter)
ax2.hist(x = data2, bins = inter)
ax3.hist(x = data3, bins = inter)
ax4.hist(x = data4, bins = inter)
ax1.set_title('RNA data set: '+ str(len(data)))
ax2.set_title('Pacientes con estudio WES: '+ str(len(data2)))
ax3.set_title('Pacientes con estudio RNA y WES: '+ str(len(data3)))
ax4.set_title('Clinic and Immune phenotype data set: '+ str(len(data4)))

plt.show()

f = open("RNA_pfs.txt", "r")
for linea in f:
   l = float(linea[:-2])
   data4.append(l)
f.close()
print('Mediana de los datos PFS (Pacientes Clinico): '+ str(statistics.median(data4)))