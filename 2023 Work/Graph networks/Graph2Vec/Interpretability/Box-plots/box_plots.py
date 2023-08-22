import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
import seaborn as sns

os.chdir('2023 Work/Graph networks/Graph2Vec/Interpretability/Box-plots/')

results = pd.read_excel('results_random_forest.xlsx')
print(results.iloc[0].values[1:5])

# Creating dataset
final_data = []; final_label = []
A = 0
N = 108
data = []
for i in range (A,N):
    data = results.iloc[i].values[1:5] # (all)
    # data.append(results.iloc[i].values[1]) # bio-niv
    # data.append(results.iloc[i].values[2]) # bio-ave-axi
    # data.append(results.iloc[i].values[3]) # ppt-niv
    # data.append(results.iloc[i].values[4]) # ppt-ave-axi
    label = results.iloc[i].values[0:1]
    final_data.append(data)
    final_label.append(label)
    data = []

fig = plt.figure(figsize =(10, 7))
 
 
# Creating plot
bp = plt.boxplot(final_data)
plt.title("ROC AUC results for RF classifier (A+A-BioGrid and A+A-PPT)")
plt.xticks(np.arange(1,N+1-A),final_label,rotation=90)
plt.subplots_adjust(bottom=0.2, left=0.02, right=0.999)
for i in np.arange(A+9.07,N,9.07):
    plt.axvline(x=i, ymin=0.0, ymax=1.0)
# show plot
X=[]
Y=[]
y=0; x = 0
for m in bp['medians']:
    [[x0, x1],[y0,y1]] = m.get_data()
    x1 = np.mean((x0,x1))
    x = x + x1
    y = y + np.mean((y0,y1))
    if x1%9==0:
        Y.append(y/9)
        X.append(x/9)
        y = 0; x = 0
plt.plot(X,Y,c='C1')
plt.xticks(fontsize=8)
plt.show()
