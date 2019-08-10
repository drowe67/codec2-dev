#!/usr/bin/python3
# phasenn_train.py
#
# David Rowe August 2019
#
# Keras model for estimating the phase of sinusoidally modelled speech

# To generate features:
#   $ ./c2sim ~/Downloads/all_speech_8k.sw --dumpphase_nnl train.f32

import numpy as np
import sys
from keras.layers import Dense
from keras import models,layers
from keras import initializers

if len(sys.argv) < 2:
    print("usage: phasenn_train.py train.f32")
    sys.exit(0)

# constants

max_amp     = 80   # sparsevector covering 0 ... Fs/2
nb_features = 243  # number of sparse features/row in input training data file
nb_epochs   = 10

# load training data

feature_file = sys.argv[1]
features = np.fromfile(feature_file, dtype='float32')
nb_frames = int(len(features)/nb_features)
print("nb_frames: %d" % (nb_frames))

# 0..80    log10(A)
# 81..161  cos(phi)
# 162..242 sin(phi)

features = np.reshape(features, (nb_frames, nb_features))
print("features shape:")
print(features.shape)

# So the idea is we can predict the next frames phases from the
# current frame, and the magnitude spectrum.  For voiced speech, the
# sinusoids are continuous, so can be predicted from frame to frame if
# you know the frequency and previous phase.  We encode the frequency
# as the position in the sprase vector.

# Cascased with that is phase spectra due to dispersion of the phase
# response of the vocal tract filter, e.g. a large dispersion around
# resonances.  We supply the magnitude spectra to help model the vocal
# tract filter phase.

# Unvoiced speech has more random phase.  Hopefully the NN can work
# out if the speech is voiced or unvoiced from the magnitide spectra.

# The phase is encoded using cos and sin of the phase, as these are
# bounded by +/-1

# So input features are this frame's log(A), and last frames phase.
# The output features we are trying to model are this frames phase.

train = np.concatenate( (features[1:,:max_amp+1], features[:-1,max_amp+1:]) )
target = features([1:,max_amp+1:])

model = models.Sequential()
model.add(layers.Dense(256, activation='relu', input_dim=nb_input))
model.add(layers.Dense(256, activation='relu'))
model.add(layers.Dense(256, activation='relu'))
model.add(layers.Dense(nb_ouput, activation='linear'))

# Custom loss function that measures difference in phase just at
# non-zero elements of target (ytrue).  This could be extended to
# weight each phase error by the (log) Amplitude of each harmonic

import keras.backend as K
def customLoss(yTrue, yPred):
    # generate a mask vector with 1's on non zero values of yTrue
    mask = abs(K.sign(yTrue))
    # collect error in cos() and sin() terms, ignoring yPred values outside of
    # harmonics we care about
    error = yTrue - mask * yPred
    return K.sum(error * error)

# Compile our model 

from keras import optimizers
model.compile(loss=customLoss, optimizer='sge')

# fit model, using 20% of our data for validation

history = model.fit(train, target, validation_split=0.2, batch_size=32, epochs=nb_epochs)
model.save("phasenn_model.h5")

import matplotlib.pyplot as plt

plot_en = 0;
if plot_en:
    plt.figure(1)
    plt.plot(10*np.sqrt(history.history['loss']))
    plt.plot(10*np.sqrt(history.history['val_loss']))
    plt.title('model loss')
    plt.ylabel('rms error (rad)')
    plt.xlabel('epoch')
    plt.legend(['train', 'valid'], loc='upper right')
    plt.show()

# run model on training data and measure variance, should be similar to training "loss"

train_out = model.predict(train)
err = (train_out - target)
var = np.var(err)
std = np.std(err)
print("var: %f std: %f" % (var,std))

