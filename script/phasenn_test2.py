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
import keras.backend as K

import matplotlib.pyplot as plt

# constants

N                 = 80      # number of time domain samples in frames
nb_samples        = 20000
nb_batch          = 512
nb_epochs         = 10
width             = 256     # frequency bins, width of input/output vectors
fo_min            = 50.0
fo_max            = 400.0
Fs                = 8000.0

# Generate training data.  Given the phase at the start of the frame,
# and the frequency, determine the phase at the end of the frame

# phase encoded as cos,sin pairs
phase_start = np.zeros((nb_samples, 2*width))
phase_end = np.zeros((nb_samples, 2*width))

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
        bin = int(Wm*width/np.pi)
        r = np.random.rand(1)
        phase_start_pol = -np.pi + r[0]*2*np.pi
        phase_start[i,2*bin] = np.cos(phase_start_pol)
        phase_start[i,2*bin+1] = np.sin(phase_start_pol)
        phase_end_pol = phase_start_pol + N*Wm
        phase_end[i,2*bin] = np.cos(phase_end_pol)
        phase_end[i,2*bin+1] = np.sin(phase_end_pol)
        
print(phase_start.shape)                       
print(phase_end.shape)

model = models.Sequential()
model.add(layers.Dense(width*20, activation='relu', input_dim=2*width))
model.add(layers.Dense(width*20, activation='relu'))
model.add(layers.Dense(2*width, activation='linear'))
model.summary()

# Compile and fit our model 

from keras import optimizers
model.compile(loss='mse', optimizer='adam')
history = model.fit(phase_start, phase_end, batch_size=nb_batch, epochs=nb_epochs)

# measure error over all samples

phase_end_est = model.predict(phase_start)
err = (phase_end_est - phase_end)
var = np.var(err)
std = np.std(err)
print("rect var: %f std: %f" % (var,std))

err_angle = np.arctan2(err[:,1::2],err[:,::2])
print(err_angle.shape)
var = np.var(err_angle)
std = np.std(err_angle)
print("angle var: %f std: %f" % (var,std))

for i in range(nb_samples):
    print("Wo: %f L: %d\n" % (Wo[i], L[i]))
    for m in range(1,L[i]+1):
        Wm = m*Wo[i]
        bin = int(np.floor(Wm*width/np.pi))
        print("bin: %3d phase_end: % 5.4f % 5.4f est: % 5.4f % 5.4f err: % 5.4f % 5.4f" %
              (bin, phase_end[i][2*bin], phase_end[i][2*bin+1],
               phase_end_est[i][2*bin], phase_end_est[i][2*bin+1], err[i][2*bin], err[i][2*bin+1]))
    if i == 1:
        exit()      

plot_en = 0;
if plot_en:
    plt.figure(1)
    plt.plot(10*np.log10(history.history['loss']))
    plt.title('model loss (dB)')
    plt.xlabel('epoch')
 
    #plt.figure(2)
    #plt.hist(err_angle, bins=20)
    #plt.show()
   
