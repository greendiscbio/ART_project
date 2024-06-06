# 1. INSTALL DEPENDENCIES
"""
# pip install ray torch torchvision
# pip install bayesian-optimization

# from google.colab import drive
# drive.mount('/content/drive')

"""
# 2. IMPORT LIBRARIES

import pandas as pd
import networkx as nx
import tensorflow as tf
import numpy as np

from sklearn.model_selection import train_test_split 
from sklearn.model_selection import StratifiedKFold

from tensorflow.keras import layers, models

from keras import layers
from keras import backend as K
from keras.callbacks import CSVLogger

from ray import tune, air
from ray.tune.schedulers.hb_bohb import HyperBandForBOHB
from ray.tune.search.bohb import TuneBOHB

# 3. INITIALIZE MAIN PARAMETERS

test_size = 0.1
kf_splits = 10
epochs = 300
samples = 200

# 4. CREATE MODEL

#def Model(config):
#  model = models.Sequential()
#  model.add(layers.Conv2D(filters=config['conv_block1_filters'], kernel_size=(5, 5), activation='relu', input_shape = (100,40,1)))
#  model.add(layers.MaxPooling2D((2, 2)))
#  model.add(layers.Conv2D(filters=config['conv_block2_filters'], kernel_size=(3, 3), activation='relu'))
#  model.add(layers.MaxPooling2D((2, 2)))
#  model.add(layers.Conv2D(filters=config['conv_block3_filters'], kernel_size=(3, 3), activation='relu'))
#  model.add(layers.MaxPooling2D((2, 2)))
#  model.add(layers.Flatten())
#  model.add(layers.Dense(units=config['fc1_units']))
#  model.add(layers.Dense(units=config['fc2_units']))
#  model.add(layers.Dense(units=config['fc3_units']))
#  model.add(layers.Dense(1, activation='sigmoid'))
#  model.summary()
#  return model


def Model(config):
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


# 5. TRAIN MODEL

def train_rna_structure(config, checkpoint_dir=None):
  data = pd.read_csv("/home/sandralonso/2023_work/spectral_clustering/spectral_clustering_feature_selection/CNN/Data/clinic_and_RNA_data_raw_NIVOLUMAB.csv")
  initial_G = nx.read_edgelist("/home/sandralonso/2023_work/spectral_clustering/spectral_clustering_feature_selection/CNN/Data/biogrid_found_genes_level_2_hvgs.edgelist")  
  G = nx.Graph(initial_G)
  print("Complete graph:")
  print(G)
  Gcc = sorted(nx.connected_components(G), key=len, reverse=True)
  G = G.subgraph(Gcc[0])
  print("Maximum component:")
  print(G)

  Y = []
  X = []
  for patient in range(len(data)): 
    print('patient: ' +str(patient))
    final_matrix = pd.read_csv('/home/sandralonso/2023_work/spectral_clustering/spectral_clustering_feature_selection/CNN/Patients_zoom/final_matrix_'+str(patient)+'.csv')
    if data.PFS[patient] < 3:
      Y.append(0)
    else:
      Y.append(1)
    if 'Unnamed: 0' in final_matrix.columns:
      final_matrix = final_matrix.drop(columns=['Unnamed: 0'])
    X.append(final_matrix)

  
  

  # Train-Test split
  X_train, X_test, Y_train, Y_test = train_test_split(X, Y, test_size=test_size, stratify=Y)
  with open("/home/sandralonso/2023_work/spectral_clustering/spectral_clustering_feature_selection/CNN/Results_structure/Test_patients_X.txt", "w") as file:
        for i in range(len(X_test)):
           file.write(str(X_test[i]) + '\n')
  with open("/home/sandralonso/2023_work/spectral_clustering/spectral_clustering_feature_selection/CNN/Results_structure/Test_patients_Y.txt", "w") as file:
        for i in range(len(Y_test)):
           file.write(str(Y_test[i]) + '\n')
  

  # Initialize result dictionary
  history_m = ({"loss":0, "val_loss": 0, "accuracy":0, "val_accuracy": 0, "f1" : 0, "val_f1": 0, "precision" : 0, "val_precision": 0, "recall" : 0, "val_recall": 0})
  
  # Cross validation
  kf=StratifiedKFold(n_splits=kf_splits, random_state=None, shuffle=False)
  current_fold = 0
  for train_index, val_index in kf.split(X_train, Y_train):
      print('Fold: '+str(current_fold))
      current_fold = current_fold + 1

      # Create  model
      model = Model(config)
      
      # Create train and validation datasets
      train_dataset = []
      val_dataset = []
      train_y_dataset = []
      val_y_dataset = []
      print("TRAIN: ", train_index, "TEST:", val_index)
      for i in train_index:
          train_dataset.append(X_train[i])
          train_y_dataset.append(Y_train[i])
      for i in val_index:
          val_dataset.append(X_train[i])
          val_y_dataset.append(Y_train[i])
      print('Patients in train dataset'+ str(len(train_dataset)))
      print('Patients in val dataset'+ str(len(val_dataset)))

      train_dataset=np.array(train_dataset)
      val_dataset=np.array(val_dataset)
      train_y_dataset=np.array(train_y_dataset)
      val_y_dataset=np.array(val_y_dataset)
            
      # Compile model with losses and metrics
      model.compile(optimizer=tf.keras.optimizers.Adam(learning_rate = config['lr']),
                  loss='binary_crossentropy',
                  metrics=['accuracy', f1, precision_m, recall_m])
  
      # Create callbacks to be used during model training
      csv_logger = CSVLogger('/home/sandralonso/2023_work/spectral_clustering/spectral_clustering_feature_selection/CNN/Results_structure/log_bo.csv', append=True, separator=';')
      
      # Start model training
      history = model.fit(train_dataset, train_y_dataset,
                          epochs=epochs,
                          batch_size = config['batch_size'],
                          validation_data=(val_dataset, val_y_dataset),
                          callbacks=[csv_logger])
      
      # Save results and add to the other folds results in order to caluclate the median value afterwards
      history_m = {
      "loss": history_m["loss"] + history.history["loss"][-1],
      "val_loss": history_m["val_loss"] + history.history["val_loss"][-1],
      "accuracy": history_m["accuracy"] + history.history["accuracy"][-1],
      "val_accuracy": history_m["val_accuracy"] + history.history["val_accuracy"][-1],
      "f1" : history_m["f1"] + history.history["f1"][-1],
      "val_f1": history_m["val_f1"] + history.history["val_f1"][-1],
      "precision" : history_m["precision"] + history.history["precision_m"][-1],
      "val_precision": history_m["val_precision"] + history.history["val_precision_m"][-1],
      "recall" : history_m["recall"] + history.history["recall_m"][-1],
      "val_recall": history_m["val_recall"] + history.history["val_recall_m"][-1]
      }

  # Calculate the median value of the results
  for i in history_m:
    history_m[i] = history_m[i] / kf_splits

  # Save the median result
  with open("/home/sandralonso/2023_work/spectral_clustering/spectral_clustering_feature_selection/CNN/Results_structure/folds_results_bo.txt", "a") as file:
    file.write(str(history_m)+'\n')

  return {"median_loss": history_m["loss"], "median_val_loss": history_m["val_loss"], "median_accuracy": history_m["accuracy"], "median_val_accuracy": history_m["val_accuracy"], "median_f1" : history_m["f1"], "median_val_f1": history_m["val_f1"], "median_precision" : history_m["precision"], "median_val_precision": history_m["val_precision"], "median_recall" : history_m["recall"], "median_val_recall": history_m["val_recall"]}

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



#search_space = {"lr": tune.choice([0.00001, 0.0001, 0.001, 0.01, 0.1]),
#                    "batch_size": tune.choice([8, 16, 32, 64, 128]), 
#                    "conv_block1_filters":tune.choice([16, 32, 64, 128, 256, 512]),
#                    "conv_block2_filters":tune.choice([16, 32, 64, 128, 256, 512]),
#                    "conv_block3_filters":tune.choice([16, 32, 64, 128, 256, 512]),
#                    "fc1_units":tune.choice([16, 32, 64, 128, 256, 512]),
#                    "fc2_units":tune.choice([16, 32, 64, 128, 256, 512]),
#                    "fc3_units":tune.choice([16, 32, 64, 128, 256, 512])
#  }


search_space = {"lr": tune.choice([0.00001, 0.0001, 0.001, 0.01, 0.1]),
                "n_fc_layers": tune.choice([0,1,2,3,4]), 
                "n_layers": tune.choice([1,2,3,4]), 
                "batch_size": tune.choice([8, 16, 32, 64, 128]), 
                "conv_block1_filters":tune.choice([16, 32, 64, 128, 256, 512]),
                "conv_block2_filters":tune.choice([16, 32, 64, 128, 256, 512]),
                "conv_block3_filters":tune.choice([16, 32, 64, 128, 256, 512]),
                "conv_block4_filters":tune.choice([16, 32, 64, 128, 256, 512]),
                "conv_block5_filters":tune.choice([16, 32, 64, 128, 256, 512]),
                "pool_type1": tune.choice(["max", "average"]),
                "pool_type2": tune.choice(["max", "average"]),
                "pool_type3": tune.choice(["max", "average"]),
                "pool_type4": tune.choice(["max", "average"]),
                "pool_type5": tune.choice(["max", "average"]),
                "fc_layer1_dropout": tune.choice([0.0, 0.1, 0.2, 0.3, 0.4, 0.5]),
                "fc_layer2_dropout": tune.choice([0.0, 0.1, 0.2, 0.3, 0.4, 0.5]),
                "fc_layer3_dropout": tune.choice([0.0, 0.1, 0.2, 0.3, 0.4, 0.5]),
                "fc_layer4_dropout": tune.choice([0.0, 0.1, 0.2, 0.3, 0.4, 0.5]),
                "fc_layer5_dropout": tune.choice([0.0, 0.1, 0.2, 0.3, 0.4, 0.5]),
                "fc_layer1_batch_norm": tune.choice([True, False]),
                "fc_layer2_batch_norm": tune.choice([True, False]),
                "fc_layer3_batch_norm": tune.choice([True, False]),
                "fc_layer4_batch_norm": tune.choice([True, False]),
                "fc_layer5_batch_norm": tune.choice([True, False]),
                "fc1_units":tune.choice([16, 32, 64, 128, 256, 512]),
                "fc2_units":tune.choice([16, 32, 64, 128, 256, 512]),
                "fc3_units":tune.choice([16, 32, 64, 128, 256, 512]),
                "fc4_units":tune.choice([16, 32, 64, 128, 256, 512]),
                "fc5_units":tune.choice([16, 32, 64, 128, 256, 512])
  }


bohb_hyperband = HyperBandForBOHB(
        time_attr="training_iteration",
        max_t=5,
    )

bohb_search = TuneBOHB()
bohb_search = tune.search.ConcurrencyLimiter(bohb_search, max_concurrent=1)

tuner = tune.Tuner(
    train_rna_structure,
    run_config=air.RunConfig(name="bohb_test", stop={"median_accuracy": 0.99}),
    tune_config=tune.TuneConfig(
        metric="median_val_f1",
        mode="max",
        scheduler=bohb_hyperband,
        search_alg=bohb_search,
        num_samples=samples,
    ),
    param_space=search_space,
)
results = tuner.fit()

print(results)
best_config = results.get_best_result().config

with open("Results_structure/best_result_bo.txt", 'w') as file:
  file.write('Best config; '+ str(best_config))
  file.write('Best result: ' +str(results.get_best_result()))
with open("Results_structure/all_results_bo.txt", 'w') as file:
  for i in results:
    file.write(str(i)+'\n')
