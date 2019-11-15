#!/usr/bin/python3
# phasenn_test9.py
#
# David Rowe Nov 2019

# Estimate an impulse position from the phase spectra of a 2nd order system excited by an impulse
#
# periodic impulse train Wo at time offset n0 -> 2nd order system -> discrete phase specta -> NN -> n0

# TODO:
#   [ ] sanity check n0 phases

import numpy as np
import sys
from keras.layers import Dense
from keras import models,layers
from keras import initializers
import matplotlib.pyplot as plt
from scipy import signal
from keras import backend as K

# constants

N                 = 80      # number of time domain samples in frame
nb_samples        = 10000
nb_batch          = 32
nb_epochs         = 25
width             = 256
pairs             = 2*width
fo_min            = 50
fo_max            = 400
Fs                = 8000

# Generate training data

amp = np.zeros((nb_samples, width))
# phase as an angle
phase = np.zeros((nb_samples, width))
# phase encoded as cos,sin pairs:
phase_rect = np.zeros((nb_samples, pairs))
Wo = np.zeros(nb_samples)
L = np.zeros(nb_samples, dtype=int)
n0 = np.zeros(nb_samples, dtype=int)
target = np.zeros((nb_samples,1))

for i in range(nb_samples):

    # distribute fo randomly on a log scale, gives us more training
    # data with low freq frames which have more harmonics and are
    # harder to match
    r = np.random.rand(1)
    log_fo = np.log10(fo_min) + (np.log10(fo_max)-np.log10(fo_min))*r[0]
    fo = 10 ** log_fo
    Wo[i] = fo*2*np.pi/Fs
    L[i] = int(np.floor(np.pi/Wo[i]))
    # pitch period in samples
    P = 2*L[i]

    r = np.random.rand(3)

    # sample 2nd order IIR filter with random peak freq (alpha) and peak amplitude (gamma)
    alpha = 0.1*np.pi + 0.8*np.pi*r[0]
    gamma = r[1]
    w,h = signal.freqz(1, [1, -2*gamma*np.cos(alpha), gamma*gamma], range(1,L[i])*Wo[i])
    
    # select n0 between 0...P-1 (it's periodic)
    n0[i] = r[2]*P
    e = np.exp(-1j*n0[i]*range(1,L[i])*Wo[i])
    
    for m in range(1,L[i]):
        bin = int(np.round(m*Wo[i]*width/np.pi))
        mWo = bin*np.pi/width
        
        amp[i,bin] = np.log10(abs(h[m-1]))
        phase[i,bin] = np.angle(h[m-1]*e[m-1])
        phase_rect[i,2*bin]   = np.cos(phase[i,bin])
        phase_rect[i,2*bin+1] = np.sin(phase[i,bin])

        # target is n0 in rec coords                      
        target[i] = n0[i]/P
                              
model = models.Sequential()
model.add(layers.Dense(pairs, activation='relu', input_dim=pairs))
model.add(layers.Dense(pairs, activation='relu'))
model.add(layers.Dense(1))
model.summary()

from keras import optimizers
sgd = optimizers.SGD(lr=0.08, decay=1e-6, momentum=0.9, nesterov=True)
model.compile(loss="mse", optimizer=sgd)
history = model.fit(phase_rect, target, batch_size=nb_batch, epochs=nb_epochs)

# measure error in rectangular coordinates over all samples

target_est = model.predict(phase_rect)
err = target - target_est
var = np.var(err)
std = np.std(err)
print("var: %f std: %f" % (var,std))
print(target.shape, target_est.shape, err.shape)

plot_en = 1;
if plot_en:
    plt.figure(1)
    plt.plot(history.history['loss'])
    plt.title('model loss')
    plt.xlabel('epoch')
    plt.show(block=False)
 
    plt.figure(2)
    plt.hist(err, bins=20)
    plt.show(block=False)

    plt.figure(3)
    plt.plot(target[:10],'b')
    plt.plot(target_est[:10],'g')
    plt.show(block=False)

    # click on last figure to close all and finish
    plt.waitforbuttonpress(0)
    plt.close()
