#!/usr/bin/python3
# phasenn_test6.py
#
# David Rowe Oct 2019

# Extending test5 to deal with input/output vectors that have slightly
# different "rates", i.e.  Wo changing across the frame which is usual
# for voiced speech.

import numpy as np
import sys
from keras.layers import Dense
from keras import models,layers
from keras import initializers
import matplotlib.pyplot as plt

# constants

N                 = 80      # number of time domain samples in frame
nb_samples        = 1000000
nb_batch          = 32
nb_epochs         = 10
width             = 256
pairs             = 2*width
fo_min            = 50
fo_max            = 400
Fs                = 8000
dfo               = 0.05

# Generate training data.  Given the phase at the start of the frame,
# and the frequency, determine the phase at the end of the frame

# phase encoded as cos,sin pairs ref:
phase_start = np.zeros((nb_samples, pairs))
phase_end = np.zeros((nb_samples, pairs))
Wo_N = np.zeros(nb_samples)
Wo_0 = np.zeros(nb_samples)
L = np.zeros(nb_samples, dtype=int)

for i in range(nb_samples):

    # parameters at time 0 (start of current frame)
    r = np.random.rand(1)
    fo_0 = fo_min + (fo_max-fo_min)*r[0]
    Wo_0[i] = fo_0*2*np.pi/Fs
    L_0 = int(np.floor(np.pi/Wo_0[i]))
 
    # parameters at time N (end of current frame), allow a df0 freq change
    # across frame, typical of voiced speech
    r = np.random.rand(1)
    fo_N = fo_0 + (-2*dfo + dfo*r[0])*fo_0
    fo_N = np.max((fo_min, fo_N))
    fo_N = np.min((fo_max, fo_N))
    #fo_N = fo_0
    Wo_N[i] = fo_N*2*np.pi/Fs
    L_N = int(np.floor(np.pi/Wo_N[i]))
    L[i] = np.min((L_0, L_N))
    #print("fo: %f %f L: %d %d min: %d" % (fo_0, fo_N, L_0, L_N, L[i]))

    for m in range(1,L[i]):
        bin_0 = int(np.round(m*Wo_0[i]*width/np.pi))
        mWo_0 = bin_0*np.pi/width
        bin_N = int(np.round(m*Wo_N[i]*width/np.pi))
        mWo_N = bin_N*np.pi/width
        #print("m: %d bin_0: %d bin_N: %d" % (m, bin_0,bin_N))
        
        r = np.random.rand(1)
        phase_start_pol = -np.pi + r[0]*2*np.pi
        phase_start[i,2*bin_0]   = np.cos(phase_start_pol)
        phase_start[i,2*bin_0+1] = np.sin(phase_start_pol)

        # phase shift average of two frequencies
        phase_end_pol = phase_start_pol + N*(mWo_0 + mWo_N)/2

        phase_end[i,2*bin_N]   = np.cos(phase_end_pol)
        phase_end[i,2*bin_N+1] = np.sin(phase_end_pol)

print(Wo_0.shape, Wo_N.shape, phase_start.shape)
input = np.column_stack([Wo_0, Wo_N, phase_start])
print(input.shape)                       
print(phase_end.shape)

model = models.Sequential()
model.add(layers.Dense(pairs, activation='relu', input_dim=(pairs+2)))
#model.add(layers.Dense(pairs, activation='relu', input_dim=pairs))
model.add(layers.Dense(pairs))
model.summary()

# Compile and fit our model 

from keras import optimizers
sgd = optimizers.SGD(lr=0.04, decay=1e-6, momentum=0.9, nesterov=True)
model.compile(loss='mse', optimizer=sgd)
history = model.fit(input, phase_end, batch_size=nb_batch, epochs=nb_epochs)

# measure error in rectangular coordinates over all samples

phase_end_est = model.predict(input)
ind = np.nonzero(phase_end)
err = (phase_end_est[ind] - phase_end[ind])
var = np.var(err)
std = np.std(err)
print("rect var: %f std: %f" % (var,std))
print("angle var: %4.2f std: %4.2f degs" % (var*180/np.pi,std*180/np.pi))
#print(err[:5,:])

# approximation of angular error (for small angles) is y coord (sin)
# or error.  We don't actually care how close the (x,y) points are,
# just the error in angle.  This implies there is probably a better
# cost function that min MSE on rect coords.

err_angle = np.arctan2(err[1::2], 1)
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

    plt.figure(3)
    plt.title('sample vectors and error')
    for r in range(12):
        plt.subplot(3,4,r+1)
        phase = np.zeros(width)
        phase_est = np.zeros(width)
        for m in range(1,L[r]):
            wm = m*Wo_N[r]
            bin = int(np.round(wm*width/np.pi))
            phase[m] = np.arctan2(phase_end[r,2*bin+1], phase_end[r,2*bin])
            phase_est[m] = np.arctan2(phase_end_est[r,2*bin+1], phase_end_est[r,2*bin])
    
        plt.plot(phase[1:L[r]+1]*180/np.pi,'g')
        #err = phase[1:L[r]+1] - phase_est[1:L[r]+1]
        #err = err + 2*np.pi*np.floor(err/(2*np.pi))
        if r == 0:
            err = phase[1:L[r]+1] - phase_est[1:L[r]+1]
            print(err)
            print(np.round(err/(2*np.pi)))
        plt.plot(phase_est[1:L[r]+1]*180/np.pi,'r')
    plt.show()
   
