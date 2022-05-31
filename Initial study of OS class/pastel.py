import numpy as np
import matplotlib.pyplot as plt

data = []
f = open("RNA_os.txt", "r")
for linea in f:
   l = float(linea[:-2])
   data.append(l)
f.close()

data4 = []
f = open("Clinico_os.txt", "r")
for linea in f:
   l = float(linea[:-2])
   data4.append(l)
f.close()

result = [len(list(filter(lambda x:(x > 20), data))), len(list(filter(lambda x:(x < 20), data)))]
result4 = [len(list(filter(lambda x:(x > 20), data4))), len(list(filter(lambda x:(x < 20), data4)))]

names = ["os > 20", "os < 20"]

fig, ((ax1, ax4)) = plt.subplots(2, 1)
fig.suptitle('Proporción PFS en meses de las distintas tablas')
ax1.pie(result, labels=names, autopct="%0.1f %%")
ax4.pie(result4, labels=names, autopct="%0.1f %%")
ax1.set_title('Pacientes con estudio RNA: '+ str(len(data)))
ax4.set_title('Pacientes solo con estudio Clinico: '+ str(len(data4)))

plt.show()