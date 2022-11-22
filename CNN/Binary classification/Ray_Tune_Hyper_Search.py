#################################### LIBRARIES##########################################

import numpy as np
import torch
import torch.optim as optim
import torch.nn as nn
import torch.nn.functional as F
import pandas as pd
import tensorflow as tf
import random
import ConfigSpace as CS

from torchvision import datasets, transforms
from torch.utils.data import DataLoader
from ray import tune
from ray.tune.schedulers import ASHAScheduler
from ray.tune.schedulers import AsyncHyperBandScheduler
from ray.tune.suggest import ConcurrencyLimiter
from keras import backend as K
from sklearn import preprocessing
from sklearn.model_selection import train_test_split
from ray.tune.schedulers.hb_bohb import HyperBandForBOHB
from ray.tune.suggest.bohb import TuneBOHB

########################################################################################


#################################### MODEL #############################################

def ConvNet(config, len_classes=2):
    input = tf.keras.layers.Input(shape=(43893, 1))
    x = input
    x = tf.keras.layers.Conv1D(filters=config['conv_block1_filters'], kernel_size=(8), strides=1)(x)
    x = tf.keras.layers.BatchNormalization()(x)
    x = tf.keras.layers.Activation('relu')(x)

    if config['fc_layer_type'] == 'dense':
        if config['pool_type'] == 'max':
            x = tf.keras.layers.GlobalMaxPooling1D()(x)
        else:
            x = tf.keras.layers.GlobalAveragePooling1D()(x)

        # Fully connected layer 1
        x = tf.keras.layers.Dense(units=config['fc1_units'])(x)
        x = tf.keras.layers.Dropout(config['dropout_rate'])(x)
        x = tf.keras.layers.BatchNormalization()(x)
        x = tf.keras.layers.Activation('relu')(x)

        # Fully connected layer 2
        x = tf.keras.layers.Dense(units=len_classes)(x)
        x = tf.keras.layers.Dropout(config['dropout_rate'])(x)
        x = tf.keras.layers.BatchNormalization()(x)
        predictions = tf.keras.layers.Dense(1,tf.keras.layers.Activation('sigmoid'))(x)

    else:
        # Fully connected layer 1
        x = tf.keras.layers.Conv1D(filters=config['fc1_units'], kernel_size=1, strides=1)(x)
        x = tf.keras.layers.Dropout(config['dropout_rate'])(x)
        x = tf.keras.layers.BatchNormalization()(x)
        x = tf.keras.layers.Activation('relu')(x)


        # Fully connected layer 2
        x = tf.keras.layers.Conv1D(filters=len_classes, kernel_size=1, strides=1)(x)
        x = tf.keras.layers.Dropout(config['dropout_rate'])(x)
        x = tf.keras.layers.BatchNormalization()(x)

        if config['pool_type'] == 'max':
            x = tf.keras.layers.GlobalMaxPooling1D()(x)
        else:
            x = tf.keras.layers.GlobalAveragePooling1D()(x)

        predictions = tf.keras.layers.Dense(1,tf.keras.layers.Activation('sigmoid'))(x)

    model = tf.keras.Model(inputs=input, outputs=predictions)

    print(model.summary())
    print(f'Total number of layers: {len(model.layers)}')

    return model

########################################################################################

#################################### METRICS ###########################################
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

def f1_m(y_true, y_pred):
    precision = precision_m(y_true, y_pred)
    recall = recall_m(y_true, y_pred)
    return 2*((precision*recall)/(precision+recall+K.epsilon()))


#################################### MODEL TRAINING######################################

def train_mnist(config):
  path ='/home/sandralonso/cnn_new/Clinical_data_and_RNA_total_Features_PFS.csv'
  data_frame = pd.read_csv(path)

  X = data_frame.iloc[:,28:43921  ]
  Y=[]
  for i in range (len(data_frame)):
      if data_frame.PFS[i]<3: # If PFS is lower than 3 months, I will consider it as NonResponder (NR)
          Y.append(0)
      else:
          Y.append(1)# If PFS is over 3 months, I will consider it as Responder (R)
  scaler = preprocessing.MinMaxScaler()
  names = X.columns
  d = scaler.fit_transform(X)
  X = pd.DataFrame(d, columns=names)
  XTrain, XTest, yTrain, yTest = train_test_split(X, Y, test_size=0.10, random_state=125, stratify = Y)
  # Convert sets to arrays
  XTrain = XTrain.values
  XTest = XTest.values
  # It is mandatory to transform Y list into array for trainning the model
  yTrain=np.array(yTrain)
  yTest=np.array(yTest)

  X_train = XTrain.reshape(XTrain.shape[0], 43893 , 1)
  X_test = XTest.reshape(XTest.shape[0], 43893, 1)
  X_train = X_train.astype('float32')
  X_test = X_test.astype('float32')
  # Create model
  model = ConvNet(config)

  # Compile model with losses and metrics
  model.compile(optimizer=tf.keras.optimizers.RMSprop(learning_rate =config['lr']),
                  loss='binary_crossentropy',
                  metrics=['accuracy', f1_m,precision_m, recall_m
                           ])

  # Start model training
  history_m = model.fit(X_train, yTrain,
                      epochs=100,
                      validation_data=(X_test, yTest))

  history_m = {
  "loss": history_m.history["loss"][0],
  "val_loss": history_m.history["val_loss"][0],
  "accuracy": history_m.history["accuracy"][0],
  "val_accuracy": history_m.history["val_accuracy"][0]
  }
  return history_m

########################################################################################

################################### CONFIG SPACE #######################################

config_space = CS.ConfigurationSpace()
config_space.add_hyperparameter(
    CS.CategoricalHyperparameter(
        name="lr", choices=[ 0.001, 0.01]))
config_space.add_hyperparameter(
    CS.CategoricalHyperparameter(
        name="dropout_rate", choices=[0.1, 0.2, 0.3, 0.4, 0.5]))
config_space.add_hyperparameter(
    CS.CategoricalHyperparameter(
        name="conv_block1_filters", choices=[8, 16, 32, 64, 128]))
config_space.add_hyperparameter(
    CS.CategoricalHyperparameter(
        name="fc_layer_type", choices= ['dense', 'convolution']))
config_space.add_hyperparameter(
    CS.CategoricalHyperparameter(
        name="pool_type", choices= ['max', 'average']))
config_space.add_hyperparameter(
    CS.CategoricalHyperparameter(
        name="fc1_units", choices=[8, 16, 32, 64, 128]))

########################################################################################

############################# RUN HYPERPARAM SEARCH ####################################

algo = TuneBOHB(
    config_space, metric="val_accuracy", mode="max", max_concurrent=1)
bohb = HyperBandForBOHB(
    time_attr="training_iteration",
    metric="val_accuracy",
    mode="max",
    max_t=2)
tune.run(train_mnist, scheduler=bohb, num_samples=100, search_alg=algo)

########################################################################################
