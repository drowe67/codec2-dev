#!/usr/bin/python3
# phasenn_test6.py
#
# David Rowe Oct 2019

# Keras model for testing phase modelling using NNs. Extending test6
# to deal with input/output vectors that have slightly different "rates", i.e.
# Wo changing across the frame which is usual for voiced speech.

# TODO: clean up nomenclature (0,N) <-> (start,end)

import numpy as np
import sys
from keras.layers import Dense
from keras import models,layers
from keras import initializers
import keras.backend as K

import matplotlib.pyplot as plt

# constants

N                 = 80      # number of time domain samples in frame
nb_samples        = 10000
nb_batch          = 1
nb_epochs         = 10
width             = 256
pairs             = 2*width
fo_min            = 50
fo_max            = 400
Fs                = 8000

# Generate training data.  Given the phase at the start of the frame,
# and the frequency, determine the phase at the end of the frame

# phase encoded as cos,sin pairs
phase_start = np.zeros((nb_samples, pairs))
phase_end = np.zeros((nb_samples, pairs))

for i in range(nb_samples):

    # parameters at time 0 (start of current frame)
    r = np.random.rand(1)
    fo_0 = fo_min + (fo_max-fo_min)*r[0]
    Wo_0 = fo_0*2*np.pi/Fs
    L_0 = int(np.floor(np.pi/Wo_0))

    # parameters at time N (end of current frame), allow a 10% freq change
    # across frame, typical of voiced speech
    r = np.random.rand(1)
    fo_N = fo_0 + 0.1*f_0*r[0]
    fo_N = np.max((fo_min, fo_N))
    fo_N = np.min((fo_max, fo_N))
    Wo_N = fo_N*2*np.pi/Fs
    L_N = int(np.floor(np.pi/Wo_N))

    for m in range(np.min(L_0,L_N)):
        bin_0 = int(np.floor(m*Wo_0*width/np.pi))
        bin_N = int(np.floor(m*Wo_N*width/np.pi))
        print("fo: %f %f bin: %d %d" % (fo_0, fo_N, bin_0, bin_N))

        r = np.random.rand(1)
        phase_start_pol = -np.pi + r[0]*2*np.pi
        phase_start[i,2*bin_0]   = np.cos(phase_start_pol)
        phase_start[i,2*bin_N+1] = np.sin(phase_start_pol)

        # phase shift average of two frequencies
        phase_end_pol = phase_start_pol + N*m*(Wo_0+Wo_N)/2

        phase_end[i,2*bin_N]   = np.cos(phase_end_pol)
        phase_end[i,2*bin_N+1] = np.sin(phase_end_pol)
        
print(phase_start.shape)                       
print(phase_end.shape)

model = models.Sequential()
model.add(layers.Dense(2*pairs, activation='relu', input_dim=pairs))
model.add(layers.Dense(2*pairs, activation='relu', input_dim=pairs))
model.add(layers.Dense(pairs, input_dim=pairs))
model.summary()

# Compile and fit our model 

from keras import optimizers
sgd = optimizers.SGD(lr=0.04, decay=1e-6, momentum=0.9, nesterov=True)
model.compile(loss='mse', optimizer=sgd)
history = model.fit(phase_start, phase_end, batch_size=nb_batch, epochs=nb_epochs)

# measure error over all samples

phase_end_est = model.predict(phase_start)
err = (phase_end_est - phase_end)
var = np.var(err)
std = np.std(err)
print("rect var: %f std: %f" % (var,std))
#print(err[:5,:])

# approximation if error is small
err_angle = np.arctan2(err[:,1], 1)
#print(err[:5,:])
print(err_angle.shape)
var = np.var(err_angle)
std = np.std(err_angle)
print("angle var: %f std: %f" % (var,std))

plot_en = 1;
if plot_en:
    plt.figure(1)
    plt.plot(history.history['loss'])
    plt.title('model loss')
    plt.xlabel('epoch')
 
    plt.figure(2)
    plt.hist(err_angle, bins=20)
    plt.title('phase angle error')
    plt.show()
   
