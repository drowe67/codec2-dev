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
nb_epochs         = 100
width             = 256
pairs             = 2*width
fo_min            = 50
fo_max            = 400
Fs                = 8000

# Generate training data.  Sparse log magnitude spectrum is input,
# phase spectrum of 2nd order system the output/target

magnitude = np.zeros((nb_samples, width))
phase = np.zeros((nb_samples, pairs))

for i in range(nb_samples):

    # choose a random fo
    r = np.random.rand(1)
    fo = fo_min + (fo_max-fo_min)*r[0]
    Wo = fo*2*np.pi/Fs
    L = int(np.floor(np.pi/Wo))

    # sample 2nd order IIR filter

    w,h = signal.freqz(1, [1,0,0.81], range(1,L)*Wo)

    # map to sparse input and ouput arrays
    
    for m in range(1,L):
        bin = int(np.floor(m*Wo*width/np.pi))
        mag[i,b] = np.log10(np.abs(h[m]))
        phase_rect = h[m]/np.abs(h[m])        
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
history = model.fit(phase_start, phase_end, batch_size=nb_batch, epochs=nb_epochs)

# measure error in rectangular coordinates over all samples

phase_est = model.predict(mag)
err = (phase_est - phase)
var = np.var(err)
std = np.std(err)
print("rect var: %f std: %f" % (var,std))
#print(err[:5,:])

# approximation of angular error (for small angles) is y coord (sin) or error.  We
# don't actually care how close the (x,y) points are, just the error in angle.  This implies
# there is probably a better cost function that min MSE on rect coords.
err_angle = np.arctan2(err[:,1], 1)
#print(err[:5,:])
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
    plt.show()
   
