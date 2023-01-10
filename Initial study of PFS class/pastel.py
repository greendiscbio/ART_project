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

result = [len(list(filter(lambda x:(x > 6), data))), len(list(filter(lambda x:(x < 6), data)))]
result2 = [len(list(filter(lambda x:(x > 6), data2))), len(list(filter(lambda x:(x < 6), data2)))]
result3 = [len(list(filter(lambda x:(x > 6), data3))), len(list(filter(lambda x:(x < 6), data3)))]
result4 = [len(list(filter(lambda x:(x > 6), data4))), len(list(filter(lambda x:(x < 6), data4)))]

names = ["PFS > 6", "PFS < 6"]

fig, ((ax1, ax4), (ax3, ax2)) = plt.subplots(2, 2)
fig.suptitle('Proporción PFS en meses de las distintas tablas')
ax1.pie(result, labels=names, autopct="%0.1f %%")
ax2.pie(result2, labels=names, autopct="%0.1f %%")
ax3.pie(result3, labels=names, autopct="%0.1f %%")
ax4.pie(result4, labels=names, autopct="%0.1f %%")
ax1.set_title('RNA data set: '+ str(len(data)))
ax2.set_title('Pacientes con estudio WES: '+ str(len(data2)))
ax3.set_title('Pacientes con estudio RNA y WES: '+ str(len(data3)))
ax4.set_title('Clinic and Immune phenotype data set: '+ str(len(data4)))


plt.show()

result = [len(list(filter(lambda x:(x > 3), data))), len(list(filter(lambda x:(x < 3), data)))]
result2 = [len(list(filter(lambda x:(x > 3), data2))), len(list(filter(lambda x:(x < 3), data2)))]
result3 = [len(list(filter(lambda x:(x > 3), data3))), len(list(filter(lambda x:(x < 3), data3)))]
result4 = [len(list(filter(lambda x:(x > 3), data4))), len(list(filter(lambda x:(x < 3), data4)))]

names = ["PFS > 3", "PFS < 3"]

fig, ((ax1, ax4), (ax3, ax2)) = plt.subplots(2, 2)
fig.suptitle('Proporción PFS en meses de las distintas tablas')
ax1.pie(result, labels=names, autopct="%0.1f %%", colors= ['crimson', 'dodgerblue'])
ax2.pie(result2, labels=names, autopct="%0.1f %%")
ax3.pie(result3, labels=names, autopct="%0.1f %%")
ax4.pie(result4, labels=names, autopct="%0.1f %%")
ax1.set_title('RNA data set: : '+ str(len(data)))
ax2.set_title('Pacientes con estudio WES: '+ str(len(data2)))
ax3.set_title('Pacientes con estudio RNA y WES: '+ str(len(data3)))
ax4.set_title('Clinic and Immune phenotype data set: '+ str(len(data4)))

plt.show()

benefit = []
f = open("benefit.txt", "r")
for linea in f:
    benefit.append(linea)
f.close()
print(benefit)

print(len(list(filter(lambda x:(x == "NCB"), benefit))))
benefit2 =[ len(list(filter(lambda x:(x == "NCB\n"), benefit))), len(list(filter(lambda x:(x == "CB\n"), benefit))), len(list(filter(lambda x:(x == "ICB\n"), benefit))),len(list(filter(lambda x:(x == "NA\n"), benefit))) ]
nombres = ["NCB","CB","ICB", "NA"]
plt.pie(benefit2, labels=nombres, autopct="%0.1f %%")
plt.axis("equal")
plt.show()
