# -*- coding: utf-8 -*-
"""Autoencoder_multiclass_200.ipynb
Original file is located at
    https://colab.research.google.com/drive/1-29FiqEn9M1V-7SwJIijfmHSIvbnLawO
"""

# !pip install Keras
# !pip install keras-layer-normalization

"""# 1. Import libraries"""

import pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split
from keras.layers import Dense
import numpy as np
np.random.seed(5)
from keras.layers import Input
from keras.models import Model
from tensorflow.keras.optimizers import SGD
from sklearn.preprocessing import MinMaxScaler
import matplotlib.pyplot as pyplot
from tensorflow.keras.layers import BatchNormalization
from keras.layers import LeakyReLU
from sklearn.metrics import precision_recall_curve
from matplotlib import pyplot as plt
import collections

from google.colab import drive
drive.mount('/content/drive')

"""# 2. Read and store data"""

data = pd.read_csv("/content/drive/MyDrive/UPM/Internship/Clinical_data_and_RNA_total_Features_PFS.csv")
data

Y = data.PFS
Y

data = data.iloc[:,28:43921] 
data

"""# 4. Data preprocessing"""

data = data.iloc[:,0:-1].apply(lambda x: (x-x.min())/ (x.max() - x.min()), axis=0)

data

label =[]
for i in Y:
  if i<3:
    label.append(0)
  elif i<6:
    label.append(1)
  else:
    label.append(2)
# X_train, X_test, Y_train, Y_test = train_test_split(data, label,test_size=0.2, random_state=125, stratify=label)

from sklearn.model_selection import StratifiedKFold # import KFold
kf=StratifiedKFold(n_splits=10, random_state=None, shuffle=False)

average = []
for train_index, test_index in kf.split(data, label):
  print("TRAIN: ", train_index, "TEST:", test_index)
  X_train, X_test = data.iloc[train_index], data.iloc[test_index]
  Y_train, Y_test =  np.array(label)[train_index.astype(int)],  np.array(label)[test_index.astype(int)]
  # np.savetxt("/content/drive/MyDrive/UPM/Internship/Autoencoder/Multiclass/43/input_data_multiclass_43.csv", X_train, delimiter=',')

  # Model -------------------------------------------
  n_inputs=X_train.shape[1]
  visible = Input(shape=(n_inputs,))
    # encoder level 1
  e = Dense(n_inputs/100)(visible)
  e = BatchNormalization()(e)
  e = LeakyReLU()(e)
    # encoder level 2
  e = Dense(200)(e)
  e = BatchNormalization()(e)
  e = LeakyReLU()(e)
    # bottleneck
  n_bottleneck = 200
  bottleneck = Dense(n_bottleneck)(e)
    # decoder, level 1
  d = Dense(200)(bottleneck)
  d = BatchNormalization()(d)
  d = LeakyReLU()(d)
    # decoder level 2
  d = Dense(n_inputs/100)(d)
  d = BatchNormalization()(d)
  d = LeakyReLU()(d)
    # output layer
  output = Dense(n_inputs, activation='linear')(d)

  model = Model(inputs=visible, outputs=output)
  # Compile the model -------------------------------------------
  model.compile(optimizer='adam', loss='mse')

  # Model training -------------------------------------------
  nits = 100 
  tam_lote = 16
  history = model.fit(X_train, X_train, epochs=nits, batch_size=tam_lote, shuffle=True, validation_data=(X_test,X_test), verbose=2)

  # Loss plots -------------------------------------------
  pyplot.plot(history.history['loss'], label='train')
  pyplot.plot(history.history['val_loss'], label='test')
  pyplot.legend()
  pyplot.show()

  # Model prediction -------------------------------------------
  X_pred = model.predict(X_test)
  MSE = np.mean(np.power(X_test-X_pred,2), axis=1)
  Y_test2=np.array(Y_test, dtype=bool)

  # Precission and recall curve plt -------------------------------------------
  precision, recall, umbral = precision_recall_curve(Y_test2, MSE)
  plt.plot(umbral, precision[1:], label="Precision",linewidth=5)
  plt.plot(umbral, recall[1:], label="Recall",linewidth=5)
  plt.title('Precision y Recall para diferentes umbrales')
  plt.xlabel('Umbral')
  plt.ylabel('Precision/Recall')
  plt.legend()
  plt.show()

  # Measuring model performance
  threshold = 0.04
  Y_pred = [1 if e > threshold else 0 for e in MSE]
  match = collections.Counter(Y_pred)[0]
  print('Final accuracy on the testing dataset: ' + str(match/len(Y_pred)))
  average.append(match/len(Y_pred))
  # np.savetxt("/content/drive/MyDrive/UPM/Internship/Autoencoder/Multiclass/43/output_data_multiclass_43.csv", X_train, delimiter=',')

print("Final average: " + str(sum(average)/len(average)))

encoder = Model(inputs = visible, outputs = bottleneck)
from keras import models    
encoder.save('/content/drive/MyDrive/UPM/Internship/Autoencoder/Multiclass/200/encoder_multiclass_200.h5')
encoder = models.load_model('/content/drive/MyDrive/UPM/Internship/Autoencoder/Multiclass/200/encoder_multiclass_200.h5')
data_encode = encoder.predict(data)
data_encode = pd.DataFrame(data_encode)
data_encode.insert(0,"Target", label)
data_encode

data_encode.to_csv("/content/drive/MyDrive/UPM/Internship/Autoencoder/Multiclass/200/encoded_data_multiclass_200.csv")