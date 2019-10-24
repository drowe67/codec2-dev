#!/usr/bin/python3
# phasenn_test7.py
#
# David Rowe Oct 2019

# Keras model for testing phase modelling using NNs. Here we try
# estimating the phase of a 2nd order system from it's amplitudes.
# This models the dispersive component of speech phase spectra, up
# until now we have been testing with the linear phase component.
# This script emulates a Hilbert Transform, note however in practice
# speech is not minimum phase so HTs have there limitations for real
# speech signals.

import numpy as np
from scipy import signal
import sys
from keras.layers import Dense
from keras import models,layers
from keras import initializers
import keras.backend as K

import matplotlib.pyplot as plt

# constants

N                 = 80      # number of time domain samples in frame
nb_samples        = 10000
nb_batch          = 64
nb_epochs         = 10
width             = 256
pairs             = 2*width
fo_min            = 50
fo_max            = 400
Fs                = 8000

# Generate training data.  Sparse log magnitude spectrum is input,
# phase spectrum of 2nd order system the output/target

mag = np.zeros((nb_samples, width))
phase = np.zeros((nb_samples, pairs))

for i in range(nb_samples):

    # choose a random fo
    r = np.random.rand(1)
    fo = fo_min + (fo_max-fo_min)*r[0]
    Wo = fo*2*np.pi/Fs
    L = int(np.floor(np.pi/Wo))

    # sample 2nd order IIR filter with random peak freq and amplitude

    r = np.random.rand(2)
    alpha = 0.1*np.pi + 0.8*np.pi*r[0]
    gamma = r[1]
    w,h = signal.freqz(1, [1, -2*gamma*np.cos(alpha), gamma*gamma], range(1,L)*Wo)

    # map to sparse input and output arrays
    
    for m in range(1,L):
        bin = int(np.floor(m*Wo*width/np.pi))
        mag[i,bin] = np.log10(np.abs(h[m-1]))
        phase_rect = h[m-1]/np.abs(h[m-1])        
        phase[i,2*bin]   = phase_rect.real
        phase[i,2*bin+1] = phase_rect.imag
        
print(mag.shape)                       
print(phase.shape)

model = models.Sequential()
model.add(layers.Dense(pairs, activation='relu', input_dim=width))
model.add(layers.Dense(pairs, input_dim=pairs))
model.summary()

# Compile and fit our model 

from keras import optimizers
sgd = optimizers.SGD(lr=0.04, decay=1e-6, momentum=0.9, nesterov=True)
model.compile(loss='mse', optimizer=sgd)
history = model.fit(mag, phase, batch_size=nb_batch, epochs=nb_epochs)

# measure error in rectangular coordinates over all samples

phase_est = model.predict(mag)
err = (phase_est - phase)
var = np.var(err)
std = np.std(err)
print("rect var: %f std: %f" % (var,std))

err_angle = np.arctan2(err[:,1], 1)
print(err_angle.shape)
var = np.var(err_angle)
std = np.std(err_angle)
print("angle var: %4.2f std: %4.2f rads" % (var,std))
print("angle var: %4.2f std: %4.2f degs" % (var*180/np.pi,std*180/np.pi))

plot_en = 1;
if plot_en:
    plt.figure(1)
    plt.plot(history.history['loss'])
    plt.title('model loss')
    plt.xlabel('epoch')
 
    plt.figure(2)
    plt.hist(err_angle*180/np.pi, bins=20)
    plt.title('phase angle error (deg)')

    fig = plt.figure(3)
    ax1 = fig.add_subplot(111)
    plt.plot(20*mag[1,:])
    ax2 = ax1.twinx()
    phase = np.unwrap(np.arctan2(phase[1::2], phase[::2]))
    plt.plot(phase, 'g')
    phase_est = np.unwrap(np.arctan2(phase_est[1::2], phase_est[::2]))
    plt.plot(phase, 'r')
    
    plt.show()
   
