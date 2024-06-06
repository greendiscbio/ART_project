# import pandas as pd
# import matplotlib.pyplot as plt
# import json
# results = pd.read_csv('2023 Work/Graph networks/Spectral_clustering/Model/Ray Tune/Data/log_bo.csv', delimiter=';')
# print(results.head())

# start = 0

# # for i in range(300):
# #     found = False
# #     print('File '+str(i))
# #     metrics = results.iloc[start:start+1750]
# #     # print(start)
# #     # print(metrics)
# #     metrics_mean = metrics.groupby(['epoch']).mean()
# #     std = metrics.groupby(['epoch'])['accuracy','f1','loss','precision_m','recall_m','val_accuracy','val_f1','val_loss','val_precision_m','val_recall_m'].std()
# #     metrics_mean.to_csv('2023 Work/Graph networks/Spectral_clustering/Model/Ray Tune/Data/epoch_'+str(int(start/1750))+'.csv')
# #     std.to_csv('2023 Work/Graph networks/Spectral_clustering/Model/Ray Tune/Data/std_'+str(int(start/1750))+'.csv')
# #     start+=1750
# #     for m in range (len(metrics_mean)-10):
# #         # if m.val_accuracy > 0.59:
# #         #     print(i)
# #         #     print(m)
# #         # print(metrics_mean['val_loss'][m], metrics_mean['val_loss'][m+5])
# #         if metrics_mean['val_loss'][m]+0.01<metrics_mean['val_loss'][m+5] and metrics_mean['val_loss'][m]+0.01<metrics_mean['val_loss'][m+10]:
# #             found = True
# #             print(m, metrics_mean['val_loss'][m], metrics_mean['val_loss'][m+5], metrics_mean['accuracy'][m])
# #             print(metrics_mean['val_f1'][m], metrics_mean['val_accuracy'][m])
# #             break
# #     if found == False:
# #         print(metrics_mean['val_f1'][174], metrics_mean['val_accuracy'][174], metrics_mean['accuracy'][174])

# metrics=[]
# # for i in [1,85,88,100,158,174,185,209,212,223,269]:
# for i in [100,158,174,185,209,212,223,269]:
#     epoch = pd.read_csv('2023 Work/Graph networks/Spectral_clustering/Model/Ray Tune/Data/epoch_'+str(i)+'.csv')
#     plt.plot(epoch['loss'], label = 'loss')
#     plt.plot(epoch['val_loss'], label = 'val_loss')
#     plt.title('Results of params: '+str(i))
#     plt.legend()
#     # plt.savefig('2023 Work/Graph networks/Spectral_clustering/Model/Ray Tune/images_filter/epoch_'+str(i)+'_aloss.png')
#     plt.clf()
#     plt.plot(epoch['f1'], label = 'f1')
#     plt.plot(epoch['val_f1'], label = 'val_f1')
#     plt.title('Results of params: '+str(i))
#     plt.legend()
#     # plt.savefig('2023 Work/Graph networks/Spectral_clustering/Model/Ray Tune/images_filter/epoch_'+str(i)+'_f1.png')
#     plt.clf()    
#     metrics.append({epoch['val_loss'][174], epoch['val_f1'][174], epoch['val_accuracy'][174]})

# print(metrics)

# # epoch = pd.read_csv('2023 Work/Graph networks/Spectral_clustering/Model/Ray Tune/Data/epoc_-1.csv')
# # plt.plot(epoch['loss'], label = 'loss')
# # plt.plot(epoch['val_loss'], label = 'val_loss')
# # plt.plot(epoch['f1'], label = 'f1')
# # plt.plot(epoch['val_f1'], label = 'val_f1')
# # plt.title('Results of params: ')
# # plt.legend()
# # plt.show()

# import os
# import json
# cont = 0
# def listdirs(folder):
#    return [d for d in os.listdir(folder) if os.path.isdir(os.path.join(folder, d))]

# path = '../../Desktop/aa/'
# for dir in listdirs(path):
#   if os.path.exists(path + '/'+str(dir)+'/progress.csv') and os.path.exists(path + '/'+str(dir)+'/result.json'):
#     # print(dir)
#     with open(path + '/' +str(dir)+'/result.json') as file:
#         data = json.load(file)
        
#         if 'median_val_f1' in data:
#             if {data['median_val_loss'], data['median_val_f1'], data['median_val_accuracy']} in metrics:
#                 print(cont)
#                 print(dir, data['median_val_loss'],data['median_val_f1'], data['median_val_accuracy'])
                
#                 cont = cont+1
# print(cont)



a = [[0.52096874], [0.5832147 ], [0.58326167], [0.60300666], [0.46125925], [0.5379347 ], [0.5488673 ], [0.5579916 ], [0.47230312], [0.5766845 ], [0.48354203], [0.4808313 ], [0.48595828], [0.57383585], [0.5627122 ], [0.55547744], [0.5646019 ], [0.47199103], [0.55028343], [0.53978807], [0.4776624 ], [0.5384696 ], [0.5359379 ], [0.51808274], [0.5027422 ], [0.4585873 ], [0.46794152], [0.45384413], [0.5504683 ], [0.46193326], [0.5100645 ], [0.5430191 ], [0.53772974], [0.56943226], [0.51756674], [0.5443729 ], [0.46428138]]
for i in a:
    if i[0] < 0.468:
        print(0)
    else:
        print(1)