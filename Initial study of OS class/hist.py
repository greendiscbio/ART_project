import numpy as np
import matplotlib.pyplot as plt
import statistics
data = []
f = open("RNA_os.txt", "r")
for linea in f:
   l = float(linea[:-2])
   data.append(l)
f.close()

print('Mediana de los datos OS (Pacientes RNA): '+ str(statistics.median(data)))

data4 = []
f = open("Clinico_os.txt", "r")
for linea in f:
   l = float(linea[:-2])
   data4.append(l)
f.close()
print('Mediana de los datos OS (Pacientes Clinico): '+ str(statistics.median(data4)))

inter = np.arange(70)

fig, (ax1, ax4) = plt.subplots(2,1)
fig.suptitle('OS en distintas tablas')
ax1.hist(x = data, bins = inter)
ax4.hist(x = data4, bins = inter)
ax1.set_title('Pacientes con estudio RNA: '+ str(len(data)))
ax4.set_title('Pacientes solo con estudio Clinico: '+ str(len(data4)))
plt.show()