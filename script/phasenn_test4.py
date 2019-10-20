#!/usr/bin/python3
# phasenn_test4.py
#
# David Rowe Oct 2019

# Keras model for testing phase modelling using NNs.  Like test3 but
# with sparse input vectors of harmonic phases, as this is how phases
# will be presented from Codec 2.

# Lessons learned:
# 1. We train slower than test3, as we have less information per samples
#    with a sparse subset of discrete harmonics
# 2. We need to quantise the frequency of each bin to get low error, ow
#    we introduce a +/- 0.5bin random error in the training data, with is
#    theoretically 8 degrees, but about 10-11 in practice.

import numpy as np
import sys
from keras.layers import Dense
from keras import models,layers
from keras import initializers
import keras.backend as K

import matplotlib.pyplot as plt

# constants

N                 = 80      # number of time domain samples in frame
nb_samples        = 100000
nb_batch          = 32
nb_epochs         = 50
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
Wo = np.zeros(nb_samples)
L = np.zeros(nb_samples, dtype=int)

# Map 0...width-1 to 0...pi
for i in range(nb_samples):
    r = np.random.rand(1)
    fo = fo_min + (fo_max-fo_min)*r[0]
    Wo[i] = fo*2*np.pi/Fs
    L[i] = int(np.floor(np.pi/Wo[i]))
    for m in range(1, L[i]):
        wm = m*Wo[i]
        bin = int(np.round(wm*width/np.pi))
        # Quantise frequency to this bin, this step is important to
        # match results of test3.  Without it we are adding +/- 0.5 bin
        # uncertainty in the frequency, which will randomise phase shifts
        # across bins
        wm = bin*np.pi/width
        r = np.random.rand(1)
        phase_start_pol = -np.pi + r[0]*2*np.pi
        phase_start[i,2*bin]   = np.cos(phase_start_pol)
        phase_start[i,2*bin+1] = np.sin(phase_start_pol)
        phase_end_pol = phase_start_pol + N*wm
        phase_end[i,2*bin]   = np.cos(phase_end_pol)
        phase_end[i,2*bin+1] = np.sin(phase_end_pol)
        
print(phase_start.shape)                       
print(phase_end.shape)

# note most of these weights whould end up being 0, as we only need
# two wieghts to map each phase rotation
model = models.Sequential()
model.add(layers.Dense(pairs, bias=False, input_dim=pairs))
model.summary()

import keras.backend as K

# Compile and fit our model 

from keras import optimizers
sgd = optimizers.SGD(lr=0.02, decay=1e-6, momentum=0.9, nesterov=True)
model.compile(loss='mse', optimizer=sgd)
history = model.fit(phase_start, phase_end, batch_size=nb_batch, epochs=nb_epochs)

# measure error over non-zero samples

phase_end_est = model.predict(phase_start)
ind = np.nonzero(phase_start)
err = (phase_end_est[ind] - phase_end[ind])
var = np.var(err)
std = np.std(err)
print("rect var: %f std: %f" % (var,std))
#print(err[:5,:])

# approximation if error is small
err_angle = np.arctan2(err[1::2], 1)
#print(err[:5,:])
print(err_angle.shape)
var = np.var(err_angle)
std = np.std(err_angle)
print("angle var: %f std: %f" % (var,std))
print("angle var: %4.2f std: %4.2f degs" % (var*180/np.pi,std*180/np.pi))

plot_en = 1;
if plot_en:
    plt.figure(1)
    plt.plot(history.history['loss'])
    plt.title('model loss')
    plt.xlabel('epoch')
 
    plt.figure(2)
    plt.hist(err_angle, bins=20)
    plt.title('phase angle error')

    plt.figure(3)
    r = 0
    phase = np.zeros(width)
    phase_est = np.zeros(width)
    for m in range(1,L[r]):
        wm = m*Wo[r]
        bin = int(np.round(wm*width/np.pi))
        #print("bin: %d %f %f" % (bin, phase_end[r,2*bin], phase_end[r,2*bin+1]))
        phase[m] = np.arctan2(phase_end[r,2*bin+1], phase_end[r,2*bin])
        phase_est[m] = np.arctan2(phase_end_est[r,2*bin+1], phase_end_est[r,2*bin])
    
    plt.plot(phase[1:L[r]+1],'g')
    plt.plot(phase_est[1:L[r]+1],'r')
    
    plt.show()

