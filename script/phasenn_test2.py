#!/usr/bin/python3
# phasenn_test2.py
#
# David Rowe Oct 2019

# Keras model for testing phase modelling using NNs.  In this model,
# we map a harmonic series of phases to a sparse frequency vector.
# Self contained, generates it's own training data.

# TODO: 1/ model error acroos freqm we know LF matter most
#       2/ do I need a custom cost function?

import numpy as np
import sys
from keras.layers import Dense
from keras import models,layers
from keras import initializers
import matplotlib.pyplot as plt

# constants

N                 = 80      # number of time domain samples in frames
nb_samples        = 100000
nb_batch          = 512
nb_epochs         = 10
nb_input          = 256     # frequency bins, width of input/output vectors
nb_output         = nb_input
fo_min            = 50.0
fo_max            = 400.0
Fs                = 8000.0

# Generate training data.  Given the phase at the start of the frame,
# and the frequency, determine the phase at the end of the frame

def quant(val, levels):
    return np.round(val*levels)/levels

phase_start = np.zeros((nb_samples, nb_input))
phase_end = np.zeros((nb_samples, nb_output))
train = np.zeros((nb_samples, 2*nb_input))
target = np.zeros((nb_samples, 2*nb_output))

Wo = np.zeros(nb_samples)
L = np.zeros(nb_samples).astype(int)
for i in range(nb_samples):
    r = np.random.rand(1)
    fo = fo_min + (fo_max-fo_min)*r[0]
    Wo[i] = 2*np.pi*fo/Fs
    L[i] = int(np.floor(np.pi/Wo[i]))
    # print("fo: %f Wo: %f L: %d\n" % (fo, Wo, L))
    # generate L random phases and build sparse vectors
    for m in range(1,L[i]+1):
        Wm = m*Wo[i]
        bin = int(np.floor(Wm*nb_input/np.pi))
        r = np.random.rand(1)
        phase_start[i][bin] = -np.pi + r[0]*2*np.pi
        #phase_start[i][bin] = 0
        phase_end[i][bin] = phase_start[i][bin] + N*Wm;
        train[i][bin] = np.cos(phase_start[i][bin])
        train[i][nb_input+bin] = np.sin(phase_start[i][bin])
        target[i][bin] = np.cos(phase_end[i][bin])
        target[i][nb_output+bin] = np.sin(phase_end[i][bin])
                     
print(train.shape)                       
print(target.shape)

model = models.Sequential()
model.add(layers.Dense(256, activation='relu', input_dim=2*nb_input))
model.add(layers.Dense(256, activation='relu', input_dim=2*nb_input))
model.add(layers.Dense(2*nb_output, activation='linear'))
model.summary()

# Compile our model 

from keras import optimizers
model.compile(loss='mse', optimizer='adam')

# fit model, using 20% of our data for validation

history = model.fit(train, target, validation_split=0.2, batch_size=nb_batch, epochs=nb_epochs)

# measure error over all samples

target_est = model.predict(train)
err = (target_est - target)
print(err.shape)
var = np.var(err)
std = np.std(err)
print("rect coord var: %f std: %f" % (var,std))

# measure error over sparse samples

train_out = model.predict(train)

err_sparse = 0
nb_sparse = 0
for i in range(nb_samples):
    print("Wo: %f L: %d\n" % (Wo[i], L[i]))
    for m in range(1,L[i]+1):
        Wm = m*Wo[i]
        bin = int(np.floor(Wm*nb_input/np.pi))
        print("bin: %3d target: % 5.4f % 5.4f est: % 5.4f % 5.4f" %
              (bin, target[i][bin], target[i][nb_input+bin], target_est[i][bin], target_est[i][nb_input+bin]))
    if i == 0:
        exit()      

plot_en = 1;
if plot_en:
    plt.figure(1)
    plt.plot(10*np.log10(history.history['loss']))
    plt.plot(10*np.log10(history.history['val_loss']))
    plt.title('model loss (dB)')
    plt.xlabel('epoch')
    plt.legend(['train', 'valid'], loc='upper right')
 
    plt.figure(2)
    plt.hist(err_angle, bins=20)
    plt.show()
   
