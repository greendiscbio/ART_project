import numpy as np
import matplotlib.pyplot as plt
a = []
b = []
c = []
f = open("Current research\Informatics paper\RNA_pfs.txt", "r")
for linea in f:
    l = float(linea[:-8])
    if l<3:
        a.append(l)
    elif l<6:
        b.append(l)
    else:
        c.append(l)
f.close()
print(a)

inter = np.arange(63)
plt.hist(x = [a,b,c], bins=inter,
color=['Blue', 'Green', 'Orange'], label=['Patients whose PFS is lower than 3 months', 'Patients whose PFS is between 3 and 6 months', 'Patients whose PFS is over 6 months'])
plt.title('RNA data set: '+ str(len(a)+len(b)+len(c)))
plt.xlabel('PFS value (months)')
plt.ylabel('Patients')
plt.legend()
plt.show()
# print('Mediana de los datos PFS (Pacientes Clinico): '+ str(statistics.median(data4)))
print(len(a)/(len(a)+len(b)+len(c)))
print((len(b)+len(c))/(len(a)+len(b)+len(c)))