#################### IMPORT LIBRARIES ####################

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

#################### SETTING PARAMETERS ####################

fold = 0
batch_size = 128
lr = 0.0001
test_size = 0.1
kf_splits = 10
epochs = 162

#################### READ DATA ####################

data = pd.read_csv("Data/clinic_and_RNA_data_raw_NIVOLUMAB.csv")

initial_G = nx.read_edgelist("Data/biogrid_found_genes_level_2_hvgs.edgelist")  
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
  final_matrix = pd.read_csv('Patients_zoom/final_matrix_'+str(patient)+'.csv')
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

def create_model():
	model = models.Sequential()
	model.add(layers.Conv2D(256, (5, 5), activation='relu', input_shape = (100,40,1)))
	model.add(layers.MaxPooling2D((2, 2)))
	model.add(layers.Conv2D(256, (3, 3), activation='relu'))
	model.add(layers.MaxPooling2D((2, 2)))
	model.add(layers.Conv2D(128, (3, 3), activation='relu'))
	model.add(layers.MaxPooling2D((2, 2)))
	model.add(layers.Flatten())
	model.add(layers.Dense(128))
	model.add(layers.Dense(32))
	model.add(layers.Dense(256))
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

for train_index, val_index in kf.split(X_train, Y_train):
    print('Fold: '+str(fold))
    train_dataset=[]
    val_dataset=[]
    train_y_dataset=[]
    val_y_dataset=[]
    print("TRAIN: ", train_index, "TEST:", val_index)
    for i in train_index:
        train_dataset.append(X_train[i])
        train_y_dataset.append(Y_train[i])
    for i in val_index:
        val_dataset.append(X_train[i])
        val_y_dataset.append(Y_train[i])
    print(len(train_dataset))
    print(len(val_dataset))

    train_y_dataset=np.array(train_y_dataset)
    val_y_dataset=np.array(val_y_dataset)
    train_dataset=np.array(train_dataset)
    val_dataset=np.array(val_dataset)
    model = create_model()
#    callbacks = [tf.keras.callbacks.EarlyStopping(monitor='accuracy', mode="max", patience=5)]
    model.compile(optimizer= tf.keras.optimizers.Adam(learning_rate= lr),
                  loss='binary_crossentropy',
                  metrics=['accuracy', f1,precision_m, recall_m])
    history = model.fit(train_dataset, train_y_dataset, epochs=epochs, validation_data=(val_dataset, val_y_dataset), batch_size = batch_size)
#, class_weight=class_weight)
#, callbacks = callbacks)
    #data = {'Train loss': history.history['loss'], 'Train Acc': history.history['accuracy'], 'Train f1-score': history.history['f1'], 'Train precission': history.history['precision_m'], 'Train recall': history.history['recall_m'],'Val loss': history.history['val_loss'], 'Val Acc': history.history['val_accuracy'], 'Val f1-score': history.history['val_f1'], 'Val precission': history.history['val_precision_m'], 'Val recall': history.history['val_recall_m'] }
    #print(data)
    #json_string = json.dumps(data)
    #with open('Data/train_results.json', 'a') as outfile:
    #    outfile.write(json_string+'\n')

#################### PLOT RESULTS ####################

    plt.plot(history.history['accuracy'])
    plt.plot(history.history['val_accuracy'])
    plt.title('model accuracy')
    plt.ylabel('accuracy')
    plt.xlabel('epoch')
    plt.legend(['train', 'val'], loc='upper left')
    plt.savefig('images/acc_train_val_fold='+str(fold)+'_batch='+str(batch_size)+'_lr='+str(lr)+'_test_size='+str(test_size)+'.png')
    plt.clf()


    plt.plot(history.history['f1'])
    plt.plot(history.history['val_f1'])
    plt.title('model f1')
    plt.ylabel('f1')
    plt.xlabel('epoch')
    plt.legend(['train', 'val'], loc='upper left')
    plt.savefig('images/f1_train_val_fold='+str(fold)+'_batch='+str(batch_size)+'_lr='+str(lr)+'_test_size='+str(test_size)+'.png')
    plt.clf()

    print(history.history['loss'][0])
    plt.plot(history.history['loss'])
    plt.plot(history.history['val_loss'])
    plt.title('model loss')
    plt.ylabel('loss')
    plt.xlabel('epoch')
    plt.legend(['train', 'val'], loc='upper left')
    plt.savefig('images/loss_train_val_fold='+str(fold)+'_batch='+str(batch_size)+'_lr='+str(lr)+'_test_size='+str(test_size)+'.png')
    plt.clf()

#################### TEST MODEL ####################

    result_eval = model.evaluate(np.asarray(X_test), np.asarray(Y_test), batch_size=batch_size)
    with open('test_results/test_results__fold='+str(fold)+'_batch='+str(batch_size)+'_lr='+str(lr)+'_test_size='+str(test_size)+'.txt', 'a') as outfile:
       for r in result_eval:
          outfile.write(str(r)+", ")
       outfile.write("\n")
    result_predict = model.predict(np.asarray(X_test))
    print(result_eval)
    print(result_predict)
    print(Y_test)
    fold = fold+1
