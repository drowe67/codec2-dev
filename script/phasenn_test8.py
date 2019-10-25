#!/usr/bin/python3
# phasenn_test8.py
#
# David Rowe Oct 2019

# Estimate phase spectra from amplitude spectra for a 2nd order IIR
# filter, just like a Hilbert Transform.

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

# Generate training data.

filter_amp = np.zeros((nb_samples, width))
# phase as an angle
filter_phase = np.zeros((nb_samples, width))
# phase encoded as cos,sin pairs:
filter_phase_rect = np.zeros((nb_samples, pairs))
Wo = np.zeros(nb_samples)
L = np.zeros(nb_samples, dtype=int)

for i in range(nb_samples):

    # distribute fo randomly on a log scale, gives us more training
    # data with low freq frames which have more harmonics and are
    # harder to match
    r = np.random.rand(1)
    log_fo = np.log10(fo_min) + (np.log10(fo_max)-np.log10(fo_min))*r[0]
    fo = 10 ** log_fo
    Wo[i] = fo*2*np.pi/Fs
    L[i] = int(np.floor(np.pi/Wo[i]))
 
    # sample 2nd order IIR filter with random peak freq

    r = np.random.rand(1)
    alpha = 0.1*np.pi + 0.8*np.pi*r[0]
    gamma = 0.8; alpha = np.pi/4
    w,h = signal.freqz(1, [1, -2*gamma*np.cos(alpha), gamma*gamma], range(1,L[i])*Wo[i])
    
    for m in range(1,L[i]):
        bin = int(np.round(m*Wo[i]*width/np.pi))
        mWo = bin*np.pi/width
        
        filter_amp[i,bin] = np.log10(abs(h[m-1]))
        filter_phase[i,bin] = np.angle(h[m-1])
        filter_phase_rect[i,2*bin]   = np.cos(filter_phase[i,bin])
        filter_phase_rect[i,2*bin+1] = np.sin(filter_phase[i,bin])

model = models.Sequential()
model.add(layers.Dense(pairs, activation='relu', input_dim=width))
model.add(layers.Dense(pairs))
model.summary()

# Compile and fit our model 

from keras import optimizers
sgd = optimizers.SGD(lr=0.04, decay=1e-6, momentum=0.9, nesterov=True)
model.compile(loss='mse', optimizer=sgd)
history = model.fit(filter_amp, filter_phase_rect, batch_size=nb_batch, epochs=nb_epochs)

# measure error in rectangular coordinates over all samples

filter_phase_rect_est = model.predict(filter_amp)
ind = np.nonzero(filter_phase_rect)
err = (filter_phase_rect[ind] - filter_phase_rect_est[ind])
var = np.var(err)
std = np.std(err)
print("rect var: %f std: %f" % (var,std))

c1 = filter_phase_rect[ind]; c1 = c1[::2] + 1j*c1[1::2]
c2 = filter_phase_rect_est[ind]; c2 = c2[::2] + 1j*c2[1::2]
err_angle = np.angle(c1 * np.conj(c2))

var = np.var(err_angle)
std = np.std(err_angle)
print("angle var: %4.2f std: %4.2f rads" % (var,std))
print("angle var: %4.2f std: %4.2f degs" % (var*180/np.pi,std*180/np.pi))

def sample_model(r):
    phase = np.zeros(width, dtype=complex)
    phase_est = np.zeros(width, dtype=complex)
    phase_err = np.zeros(width, dtype=complex)
    phase_filt = np.zeros(width)
    amp_filt = np.zeros(width)
    
    for m in range(1,L[r]):
        wm = m*Wo[r]
        bin = int(np.round(wm*width/np.pi))
        phase[m] = filter_phase_rect[r,2*bin] + 1j*filter_phase_rect[r,2*bin+1]
        phase_est[m] = filter_phase_rect_est[r,2*bin] + 1j*filter_phase_rect_est[r,2*bin+1]
        phase_err[m] = phase[m] * np.conj(phase_est[m])
        amp_filt[m] = filter_amp[r,bin]
    return phase, phase_err, amp_filt
    
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
    plt.hist(Wo*(Fs/2)/np.pi, bins=20)
    plt.title('phase angle error (deg) and fo (Hz)')
    plt.show(block=False)

    plt.figure(3)
    plt.title('sample vectors and error')
    for r in range(12):
        plt.subplot(3,4,r+1)
        phase, phase_err, amp_filt = sample_model(r)    
        plt.plot(np.angle(phase[1:L[r]])*180/np.pi,'g')
        plt.plot(np.angle(phase_err[1:L[r]])*180/np.pi,'r')
        plt.ylim(-180,180)
    plt.show(block=False)

    plt.figure(4)
    plt.title('filter amplitudes')
    for r in range(12):
        plt.subplot(3,4,r+1)
        phase, phase_err, amp_filt = sample_model(r)    
        plt.plot(amp_filt[1:L[r]],'g')
    plt.show(block=False)
    
    # click on last figure to close all and finish
    plt.waitforbuttonpress(0)
    plt.close()
