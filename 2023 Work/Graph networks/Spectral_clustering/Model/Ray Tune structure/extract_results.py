import pandas as pd
import matplotlib.pyplot as plt
import json
results = pd.read_csv('2023 Work/Graph networks/Spectral_clustering/Model/Ray Tune structure/Data/log_bo.csv', delimiter=';')
print(results.head())

start = 0
metrics_=[]
for i in [4,10,11,14,15,22,23,34]:
    found = False
    print('File '+str(i))
    metrics = results.iloc[start:start+3000]
    # print(start)
    # print(metrics)
    metrics_mean = metrics.groupby(['epoch']).mean()
    std = metrics.groupby(['epoch'])['accuracy','f1','loss','precision_m','recall_m','val_accuracy','val_f1','val_loss','val_precision_m','val_recall_m'].std()
    metrics_mean.to_csv('2023 Work/Graph networks/Spectral_clustering/Model/Ray Tune structure/Results/epoch_'+str(int(start/3000))+'.csv')
    std.to_csv('2023 Work/Graph networks/Spectral_clustering/Model/Ray Tune structure/Results/std_'+str(int(start/3000))+'.csv')
    start+=3000
    
    epoch = pd.read_csv('2023 Work/Graph networks/Spectral_clustering/Model/Ray Tune structure/Results/epoch_'+str(i)+'.csv')
    plt.plot(epoch['loss'], label = 'loss')
    plt.plot(epoch['val_loss'], label = 'val_loss')
    plt.title('Results of params: '+str(i))
    plt.legend()
    plt.savefig('2023 Work/Graph networks/Spectral_clustering/Model/Ray Tune structure/Results/images_filtered/epoch_'+str(i)+'_aloss.png')
    plt.clf()
    plt.plot(epoch['f1'], label = 'f1')
    plt.plot(epoch['val_f1'], label = 'val_f1')
    plt.title('Results of params: '+str(i))
    plt.legend()
    plt.savefig('2023 Work/Graph networks/Spectral_clustering/Model/Ray Tune structure/Results/images_filtered/epoch_'+str(i)+'_f1.png')
    plt.clf()    
    plt.plot(epoch['accuracy'], label = 'f1')
    plt.plot(epoch['val_accuracy'], label = 'val_f1')
    plt.title('Results of params: '+str(i))
    plt.legend()
    plt.savefig('2023 Work/Graph networks/Spectral_clustering/Model/Ray Tune structure/Results/images_filtered/epoch_'+str(i)+'_accuracy.png')
    plt.clf() 
    metrics_.append( epoch['val_accuracy'][299])

print(metrics_)
    # for m in range (len(metrics_mean)-10):
    #     # if m.val_accuracy > 0.59:
    #     #     print(i)
    #     #     print(m)
    #     # print(metrics_mean['val_loss'][m], metrics_mean['val_loss'][m+5])
    #     if metrics_mean['val_loss'][m]+0.01<metrics_mean['val_loss'][m+5] and metrics_mean['val_loss'][m]+0.01<metrics_mean['val_loss'][m+10]:
    #         found = True
    #         print(m, metrics_mean['val_loss'][m], metrics_mean['val_loss'][m+5], metrics_mean['accuracy'][m])
    #         print(metrics_mean['val_f1'][m], metrics_mean['val_accuracy'][m])
    #         break
    # if found == False:
    #     print(metrics_mean['val_f1'][174], metrics_mean['val_accuracy'][174], metrics_mean['accuracy'][174])