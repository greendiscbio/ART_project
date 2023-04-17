#################### IMPORT LIBRARIES ####################
from sklearn.metrics import roc_auc_score
import pandas as pd
import networkx as nx
from sklearn.model_selection import train_test_split 
import tensorflow as tf
import numpy as np
from keras import layers
from keras import backend as K
import tensorflow as tf
from tensorflow.keras import layers, models
from sklearn.model_selection import StratifiedKFold
from matplotlib import pyplot as plt
import torch
from keras import layers
from keras import backend as K
from matplotlib import pyplot as plt
from sklearn import metrics
import math
#################### SETTING PARAMETERS ####################
config ={
  "batch_size": 16,
  "conv_block1_filters": 256,
  "conv_block2_filters": 256,
  "conv_block3_filters": 16,
  "conv_block4_filters": 128,
  "conv_block5_filters": 64,
  "fc1_units": 128,
  "fc2_units": 32,
  "fc3_units": 256,
  "fc4_units": 128,
  "fc5_units": 64,
  "fc_layer1_batch_norm": True,
  "fc_layer1_dropout": 0.1,
  "fc_layer2_batch_norm": False,
  "fc_layer2_dropout": 0.2,
  "fc_layer3_batch_norm": False,
  "fc_layer3_dropout": 0.2,
  "fc_layer4_batch_norm": False,
  "fc_layer4_dropout": 0.3,
  "fc_layer5_batch_norm": False,
  "fc_layer5_dropout": 0.0,
  "lr": 0.001,
  "n_fc_layers": 2,
  "n_layers": 1,
  "pool_type1": "average",
  "pool_type2": "max",
  "pool_type3": "max",
  "pool_type4": "max",
  "pool_type5": "max"
}


fold = 0
test_size = 0.2
kf_splits = 10
epochs = 450

#################### READ DATA ####################

data = pd.read_csv("../Data/clinic_and_RNA_data_raw_NIVOLUMAB.csv")

initial_G = nx.read_edgelist("../Data/biogrid_found_genes_level_2_hvgs.edgelist")  
G = nx.Graph(initial_G)
print("Complete graph description:")
print(G)
Gcc = sorted(nx.connected_components(G), key=len, reverse=True)
G = G.subgraph(Gcc[0])
print("Maximum component description:")
print(G)

Y = []
X = []
for patient in range(len(data)): # Recorro pacientes
  print('patient: ' +str(patient))
  final_matrix = pd.read_csv('../Patients_zoom/final_matrix_'+str(patient)+'.csv')
  if data.PFS[patient] < 3:
    Y.append(0)
  else:
    Y.append(1)
  if 'Unnamed: 0' in final_matrix.columns:
    final_matrix = final_matrix.drop(columns=['Unnamed: 0'])
  X.append(final_matrix)
print(X[0])
print(Y)

#################### DATA SPLIT AND STRATIFICATION ####################

X_train, X_test, Y_train, Y_test = train_test_split(X, Y, test_size=test_size, stratify=Y)
kf=StratifiedKFold(n_splits=kf_splits, random_state=None, shuffle=False)

#################### CREATE MODEL ####################

# def create_model():
# 	model = models.Sequential()
# 	model.add(layers.Conv2D(conv_block1_filters, (5, 5), activation='relu', input_shape = (100,40,1)))
# 	model.add(layers.MaxPooling2D((2, 2)))
# 	model.add(layers.Conv2D(conv_block2_filters, (3, 3), activation='relu'))
# 	model.add(layers.MaxPooling2D((2, 2)))
# 	#model.add(layers.Conv2D(conv_block3_filters, (3, 3), activation='relu'))
# 	#model.add(layers.MaxPooling2D((2, 2)))
# 	model.add(layers.Flatten())
# 	#model.add(layers.Dense(fc1_units))
# 	#model.add(layers.Dense(fc2_units))
# 	model.add(layers.Dense(fc3_units))
# 	model.add(layers.Dense(1, activation='sigmoid'))
# 	model.summary()
# 	return model

def create_model(config):
  n_layer = 1
  n_fc_layer = 0
  model = models.Sequential()
  model.add(layers.Conv2D(filters=config['conv_block1_filters'], kernel_size=(5, 5), activation='relu', input_shape = (100,40,1)))

  # Convolutional layers
  while n_layer < config['n_layers']:
    n_layer = n_layer + 1
    model.add(layers.Conv2D(filters=config['conv_block'+str(n_layer)+'_filters'], kernel_size=(3, 3), activation='relu'))
    if config['pool_type'+str(n_layer)] == 'max':
      model.add(layers.MaxPooling2D((2, 2)))
    else: 
      model.add(layers.AveragePooling2D((2, 2)))
  model.add(layers.Flatten())
  
  # Fully connected layers
  while n_fc_layer < config['n_fc_layers']:
    n_fc_layer = n_fc_layer + 1
    model.add(layers.Dense(config['fc'+str(n_fc_layer)+'_units']))
    if config['fc_layer'+str(n_fc_layer)+'_dropout'] > 0.0:
      model.add(layers.Dropout(config['fc_layer'+str(n_fc_layer)+'_dropout']))
    if config['fc_layer'+str(n_fc_layer)+'_batch_norm'] == True:
      model.add(layers.BatchNormalization())
    
  model.add(layers.Dense(1, activation='sigmoid'))

  model.summary()
  return model

#################### CREATE METRICS ####################

def recall_m(y_true, y_pred):
    true_positives = K.sum(K.round(K.clip(y_true * y_pred, 0, 1)))
    possible_positives = K.sum(K.round(K.clip(y_true, 0, 1)))
    recall = true_positives / (possible_positives + K.epsilon())
    return recall

def precision_m(y_true, y_pred):
    true_positives = K.sum(K.round(K.clip(y_true * y_pred, 0, 1)))
    predicted_positives = K.sum(K.round(K.clip(y_pred, 0, 1)))
    precision = true_positives / (predicted_positives + K.epsilon())
    return precision

def f1(y_true, y_pred):
    precision = precision_m(y_true, y_pred)
    recall = recall_m(y_true, y_pred)
    return 2*((precision*recall)/(precision+recall+K.epsilon()))

#################### TRAIN MODEL ####################

X_train=np.array(X_train)
Y_train=np.array(Y_train)
model = create_model()
model.compile(optimizer= tf.keras.optimizers.Adam(learning_rate= lr),
                loss='binary_crossentropy',
                metrics=['accuracy', f1, precision_m, recall_m])
history = model.fit(X_train, Y_train, epochs=epochs, batch_size = batch_size)


#################### PLOT RESULTS ####################

plt.plot(history.history['accuracy'])
#plt.plot(history.history['val_accuracy'])
plt.title('model accuracy')
plt.ylabel('accuracy')
plt.xlabel('epoch')
plt.legend(['train'], loc='upper left')
plt.savefig('images/acc_train_val_fold='+str(fold)+'_batch='+str(batch_size)+'_lr='+str(lr)+'_test_size='+str(test_size)+'.png')
plt.clf()


plt.plot(history.history['f1'])
#plt.plot(history.history['val_f1'])
plt.title('model f1')
plt.ylabel('f1')
plt.xlabel('epoch')
plt.legend(['train'], loc='upper left')
plt.savefig('images/f1_train_val_fold='+str(fold)+'_batch='+str(batch_size)+'_lr='+str(lr)+'_test_size='+str(test_size)+'.png')
plt.clf()

print(history.history['loss'][0])
plt.plot(history.history['loss'])
#plt.plot(history.history['val_loss'])
plt.title('model loss')
plt.ylabel('loss')
plt.xlabel('epoch')
plt.legend(['train'], loc='upper left')
plt.savefig('images/loss_train_val_fold='+str(fold)+'_batch='+str(batch_size)+'_lr='+str(lr)+'_test_size='+str(test_size)+'.png')
plt.clf()

#################### TEST MODEL ####################

result_eval = model.evaluate(np.asarray(X_test), np.asarray(Y_test), batch_size=batch_size)
with open('test_results/test_results_fold='+str(fold)+'_batch='+str(batch_size)+'_lr='+str(lr)+'_test_size='+str(test_size)+'.txt', 'a') as outfile:
   for r in result_eval:
       outfile.write(str(r)+", ")
   outfile.write("\n")
result_predict = model.predict(np.asarray(X_test))
print(result_eval)
print(result_predict)
print(Y_test)
fold = fold+1

auc = roc_auc_score(Y_test, result_predict)
print(auc)
fpr, tpr, thresholds = metrics.roc_curve(Y_test, result_predict)
#print(fpr,tpr,thresholds)
plt.plot(fpr, tpr, marker='.', label='Model')
plt.plot([0,1], [0,1], linestyle='--', label='No Skill')
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.legend()
plt.savefig('images/roc_curve.png')
plt.clf()

J = tpr - fpr
ix = np.argmax(J)
best_thresh = thresholds[ix]
print('Best Threshold=%f' % (best_thresh))
plt.plot(fpr, tpr, marker='.', label='Model')
plt.plot([0,1], [0,1], linestyle='--', label='No Skill')
plt.scatter(fpr[ix], tpr[ix], marker='o', color='black', label='Best')
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.legend()
plt.savefig('images/roc_curve_after.png')
plt.clf()


thresholds = np.arange(0, 1, 0.001)
# apply threshold to positive probabilities to create labels
def to_labels(pos_probs, threshold):
 return (pos_probs >= threshold).astype('int')


# evaluate each threshold
scores = [metrics.f1_score(Y_test, to_labels(result_predict, t)) for t in thresholds]
scores2 = [metrics.accuracy_score(Y_test, to_labels(result_predict, t)) for t in thresholds]
# get best threshold
ix = np.argmax(scores)
ix2 = np.argmax(scores2)
print('Threshold=%.3f, F-Score=%.5f' % (thresholds[ix], scores[ix]))
print('Threshold=%.3f, Accuracy=%.5f' % (thresholds[ix2], scores2[ix2]))
