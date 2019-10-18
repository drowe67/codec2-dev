#!/usr/bin/python3
# phasenn_test4.py
#
# David Rowe Oct 2019

# Keras model for testing phase modelling using NNs.  Like test3 but
# with sparse input vectors of harmonic phases, as this is how phases
# will be presented from Codec 2.

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

# Map 0...width-1 to 0...pi
for i in range(nb_samples):
    r = np.random.rand(1)
    fo = fo_min + (fo_max-fo_min)*r[0]
    Wo = fo*2*np.pi/Fs
    L = int(np.floor(np.pi/Wo))
    for m in range(L):
        wm = m*Wo
        bin = int(np.floor(wm*width/np.pi))
        r = np.random.rand(1)
        phase_start_pol = -np.pi + r[0]*2*np.pi
        phase_start[i,2*bin]   = np.cos(phase_start_pol)
        phase_start[i,2*bin+1] = np.sin(phase_start_pol)
        phase_end_pol = phase_start_pol + N*wm
        phase_end[i,2*bin]   = np.cos(phase_end_pol)
        phase_end[i,2*bin+1] = np.sin(phase_end_pol)
        
print(phase_start.shape)                       
print(phase_end.shape)

# note most of these weights whould end up being 0, as we only need two wieghts to map each phase
# rotation
model = models.Sequential()
model.add(layers.Dense(pairs, bias=False, input_dim=pairs))
model.summary()

# Compile and fit our model 

from keras import optimizers
sgd = optimizers.SGD(lr=0.02, decay=1e-6, momentum=0.9, nesterov=True)
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
   
