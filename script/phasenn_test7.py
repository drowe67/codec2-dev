#!/usr/bin/python3
# phasenn_test7.py
#
# David Rowe Oct 2019

# Extending PhaseNN model to have a linear and dispersive phase
# component.  The dispersive component is from a simple 2nd order IIR
# filter.

# TODO: [ ] remove df0 support
#       [ ] remove Wo_0, Wo_N
#       [ ] replace random starting phases with phase realted to t0

import numpy as np
import sys
from keras.layers import Dense
from keras import models,layers
from keras import initializers
import matplotlib.pyplot as plt
from scipy import signal

# constants

N                 = 80      # number of time domain samples in frame
nb_samples        = 100000
nb_batch          = 32
nb_epochs         = 25
width             = 256
pairs             = 2*width
fo_min            = 50
fo_max            = 400
Fs                = 8000
dfo               = 0.00   # dfo models aren't solved yet

# Generate training data.  Given the phase at the start of the frame,
# and the frequency, determine the phase at the end of the frame

# phase encoded as cos,sin pairs ref:
phase_start = np.zeros((nb_samples, pairs))
phase_end = np.zeros((nb_samples, pairs))
filter_amp = np.zeros((nb_samples, width))
filter_phase = np.zeros((nb_samples, width))
Wo_N = np.zeros(nb_samples)
Wo_0 = np.zeros(nb_samples)
L = np.zeros(nb_samples, dtype=int)

for i in range(nb_samples):

    # parameters at time 0 (start of current frame)
    # distribute fo randomnly on a log scale
    r = np.random.rand(1)
    log_fo_0 = np.log10(fo_min) + (np.log10(fo_max)-np.log10(fo_min))*r[0]
    fo_0 = 10 ** log_fo_0
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

    # sample 2nd order IIR filter with random peak freq and amplitude

    r = np.random.rand(2)
    alpha = 0.1*np.pi + 0.8*np.pi*r[0]
    gamma = 0.8
    #alpha = np.pi/4; gamma = 0.8
    w,h = signal.freqz(1, [1, -2*gamma*np.cos(alpha), gamma*gamma], range(1,L[i])*Wo_0[i])
    
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
        phase_end_pol = phase_start_pol + N*(mWo_0 + mWo_N)/2 + np.angle(h[m-1])

        phase_end[i,2*bin_N]   = np.cos(phase_end_pol)
        phase_end[i,2*bin_N+1] = np.sin(phase_end_pol)

        filter_amp[i,bin_0] = np.log10(abs(h[m-1]))
        filter_phase[i,bin_0] = np.angle(h[m-1])
        
input = np.column_stack([filter_amp, phase_start])

model = models.Sequential()
model.add(layers.Dense(pairs, activation='relu', input_dim=(width+pairs)))
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
err = (phase_end[ind] - phase_end_est[ind])
var = np.var(err)
std = np.std(err)
print("rect var: %f std: %f" % (var,std))

print(phase_end_est.shape, err.shape)
c1 = phase_end[ind]; c1 = c1[::2] + 1j*c1[1::2]
c2 = phase_end_est[ind]; c2 = c2[::2] + 1j*c2[1::2]
err_angle = np.angle(c1 * np.conj(c2))

print(err_angle[:5],err_angle.shape)

var = np.var(err_angle)
std = np.std(err_angle)
print("angle var: %4.2f std: %4.2f rads" % (var,std))
print("angle var: %4.2f std: %4.2f degs" % (var*180/np.pi,std*180/np.pi))

def sample_model(r):
    phase = np.zeros(width, dtype=complex)
    phase_est = np.zeros(width, dtype=complex)
    phase_err = np.zeros(width, dtype=complex)
    phase_filt = np.zeros(width)
    
    for m in range(1,L[r]):
        wm = m*Wo_N[r]
        bin = int(np.round(wm*width/np.pi))
        phase[m] = phase_end[r,2*bin] + 1j*phase_end[r,2*bin+1]
        phase_est[m] = phase_end_est[r,2*bin] + 1j*phase_end_est[r,2*bin+1]
        phase_err[m] = phase[m] * np.conj(phase_est[m])
        phase_filt[m] = filter_phase[r,bin]
    return phase, phase_err, phase_filt
    
plot_en = 1;
if plot_en:
    plt.figure(1)
    plt.plot(history.history['loss'])
    plt.title('model loss')
    plt.xlabel('epoch')
    plt.show(block=False)
 
    plt.figure(2)
    plt.subplot(211)
    plt.hist(err_angle*180/np.pi, bins=20)
    plt.subplot(212)
    plt.hist(Wo_0*(Fs/2)/np.pi, bins=20)
    plt.title('phase angle error (deg) and fo (Hz)')
    plt.show(block=False)

    plt.figure(3)
    plt.title('sample vectors and error')
    for r in range(12):
        plt.subplot(3,4,r+1)
        phase, phase_err, phase_filt = sample_model(r)    
        plt.plot(np.angle(phase[1:L[r]+1])*180/np.pi,'g')
        plt.plot(np.angle(phase_err[1:L[r]+1])*180/np.pi,'r')
        plt.plot(phase_filt[1:L[r]+1]*180/np.pi,'b')
        plt.ylim(-180,180)
    plt.show(block=False)

    # click on last figure to close all and finish
    plt.waitforbuttonpress(0)
    plt.close()
