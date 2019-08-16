#!/usr/bin/python3
# phasenn_test1.py
#
# David Rowe August 2019

# Keras model for testing phase modelling using NNs.  Self contained, generates
# it's own training data

import numpy as np
import sys
from keras.layers import Dense
from keras import models,layers
from keras import initializers

# constants

N                 = 80    # number of samples in frames
nb_samples        = 1000
nb_epochs         = 10
np_input          = 3
np_output         = 2

# Generate training data.  Given the phase at the start of the frame,
# and the frequency, determine the phase at the end of the frame

phase_start = np.zeros(nb_samples)
phase_end = np.zeros(nb_samples)
for i in range(nb_samples):
    x=np.random.rand(2)
    w = x[0]*np.pi
    phase_start[i] = x[1]*np.pi
    phase_end[i] = phase_start[i] + N*w
                     
train  = np.concatenate((w, cos(phase_start), sin(phase_start))
target = np.concatenate(cos(phase_end), sin(phase_end))
print(train.shape)                       
print(target.shape)
                        
model = models.Sequential()
model.add(layers.Dense(256, activation='relu', input_dim=np_input))
model.add(layers.Dense(256, activation='relu'))
model.add(layers.Dense(nb_output, activation='linear'))

# Compile our model 

from keras import optimizers
model.compile(loss='mse', optimizer='adam')

# fit model, using 20% of our data for validation

history = model.fit(train, target, validation_split=0.2, batch_size=32, epochs=nb_epochs)

import matplotlib.pyplot as plt

plot_en = 1;
if plot_en:
    plt.figure(1)
    plt.plot(np.sqrt(history.history['loss']))
    plt.plot(np.sqrt(history.history['val_loss']))
    plt.title('model loss')
    plt.ylabel('rms error')
    plt.xlabel('epoch')
    plt.legend(['train', 'valid'], loc='upper right')
    plt.show()

