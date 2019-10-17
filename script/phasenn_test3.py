#!/usr/bin/python3
# phasenn_test4.py
#
# David Rowe Oct 2019

# Keras model for testing phase modelling using NNs.  Extending test3
# to sparse model, which is now codec 2 phases are generated

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
nb_epochs         = 5
width             = 256
pairs             = 2*width

# Generate training data.  Given the phase at the start of the frame,
# and the frequency, determine the phase at the end of the frame

# phase encoded as cos,sin pairs
phase_start = np.zeros((nb_samples, pairs))
phase_end = np.zeros((nb_samples, pairs))

# Map 0...width-1 to 0...pi
w = np.pi/width;
for i in range(nb_samples):
    for m in range(width):
        r = np.random.rand(1)
        phase_start_pol = -np.pi + r[0]*2*np.pi
        phase_start[i,2*m]   = np.cos(phase_start_pol)
        phase_start[i,2*m+1] = np.sin(phase_start_pol)
        phase_end_pol = phase_start_pol + N*w*m
        phase_end[i,2*m]   = np.cos(phase_end_pol)
        phase_end[i,2*m+1] = np.sin(phase_end_pol)
        
print(phase_start.shape)                       
print(phase_end.shape)

# note most of these weights whould end up being 0, as we only need two wieghts to map each phase
# rotation
model = models.Sequential()
model.add(layers.Dense(pairs, bias=False, input_dim=pairs))
model.summary()

# Compile and fit our model 

from keras import optimizers
sgd = optimizers.SGD(lr=0.01, decay=1e-6, momentum=0.9, nesterov=True)
model.compile(loss='mse', optimizer=sgd)
history = model.fit(phase_start, phase_end, batch_size=nb_batch, epochs=nb_epochs)

# measure error over all samples

phase_end_est = model.predict(phase_start)
err = (phase_end_est - phase_end)
var = np.var(err)
std = np.std(err)
print("rect var: %f std: %f" % (var,std))
#print(err[:5,:])

err_angle = np.arctan2(phase_end_est[:,1], phase_end_est[:,0]) - np.arctan2(phase_end[:,1], phase_end[:,0])
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
   
