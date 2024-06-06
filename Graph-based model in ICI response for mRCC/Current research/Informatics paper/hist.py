import numpy as np
import matplotlib.pyplot as plt
a = []
b = []
c = []
all =[]
f = open("Current research\Informatics paper\RNA_pfs.txt", "r")
for linea in f:
    l = float(linea[:-8])
    all.append(l)
    if l<3:
        a.append(l)
    elif l<6:
        b.append(l)
    else:
        c.append(l)
f.close()
print(a)

plt.hist(x = [a,b,c], bins=7,

color=['Blue', 'Green', 'Orange'], label=['Patients whose PFS is lower than 3 months (NR)', 'Patients whose PFS is between 3 and 6 months (IR)', 'Patients whose PFS is over 6 months (R)'])
plt.title('RNA data: '+ str(len(a)+len(b)+len(c))+ ' patients')
plt.xlabel('PFS value (months)')
plt.ylabel('Patients')
plt.legend()
plt.show()
# print('Mediana de los datos PFS (Pacientes Clinico): '+ str(statistics.median(data4)))
print('Percentage of class 0: '+str(len(a)/(len(a)+len(b)+len(c))))
print('Percentage of class 1: '+str(len(b)/(len(a)+len(b)+len(c))))
print('Percentage of class 2: '+str(len(c)/(len(a)+len(b)+len(c))))
print('Percentage of class 1 and 2: '+str((len(b)+len(c))/(len(a)+len(b)+len(c))))

print('Patients of class 0: '+str(len(a)))
print('Patients of class 1: '+str(len(b)))
print('Patients of class 2: '+str(len(c)))
print(np.mean(all))