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
import matplotlib.pyplot as plt

# constants

N                 = 80      # number of samples in frames
nb_samples        = 100000
nb_epochs         = 20
np_input          = 3
nb_output         = 2

# Generate training data.  Given the phase at the start of the frame,
# and the frequency, determine the phase at the end of the frame

def quant(val, levels):
    return np.round(val*levels)/levels

w  = np.zeros(nb_samples)
phase_start = np.zeros(nb_samples)
phase_end = np.zeros(nb_samples)
for i in range(nb_samples):
    x=np.random.rand(2)
    #w[i] = np.pi/4 + 0.1*x[0]*np.pi
    w[i] = np.pi/8 + x[0]*np.pi/8
    phase_start[i] = -np.pi + x[1]*2*np.pi
    phase_end[i] = phase_start[i] + N*w[i];
                     
train  = np.column_stack((w, np.cos(phase_start), np.sin(phase_start)))
target = np.column_stack((np.cos(phase_end), np.sin(phase_end)))
print(train.shape)                       
print(target.shape)
                        
model = models.Sequential()
model.add(layers.Dense(256, activation='relu', input_dim=np_input))
model.add(layers.Dense(256, activation='relu', input_dim=np_input))
model.add(layers.Dense(nb_output, activation='linear'))
model.summary()

# Compile our model 

from keras import optimizers
model.compile(loss='mse', optimizer='adam')

# fit model, using 20% of our data for validation

history = model.fit(train, target, validation_split=0.2, batch_size=32, epochs=nb_epochs)

# test actual error in angle over training data

train_out = model.predict(train)
err = (train_out - target)
var = np.var(err)
std = np.std(err)
print("rect coord var: %f std: %f" % (var,std))

err_angle = np.arctan2(np.sin(err[:,0]),np.cos(err[:,1]))
var = np.var(err_angle)
std = np.std(err_angle)
print("angle var: %f std: %f" % (var,std))

plot_en = 1;
if plot_en:
    plt.figure(1)
    plt.plot(10*np.log10(history.history['loss']))
    plt.plot(10*np.log10(history.history['val_loss']))
    plt.title('model loss (dB)')
    plt.xlabel('epoch')
    plt.legend(['train', 'valid'], loc='upper right')
 
    plt.figure(2)
    plt.hist(err_angle)
    plt.show()
   
    #plt.figure(2)
    #plt.hist(target)
    #plt.show()
