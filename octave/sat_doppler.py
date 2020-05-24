#!/usr/bin/env python
#
#   Generate doppler data to be applied to a test modem waveform.
#
#   $ python sat_doppler.py
#   This will generate doppler.npz, which will be used in the fsk_demod_BER_test.py script
#   if DOPPLER_ENABLED is True.
#
#   Dependencies:
#       - python-dateutil
#       - pyephem
#       - numpy
#       - matplotlib
#
#   Copyright (C) 2020  Mark Jessop <vk5qi@rfhead.net>
#   Released under GNU GPL v3 or later

import datetime
import json
import logging
import os
import time
import traceback
import subprocess
import sys

import ephem
import numpy as np
import matplotlib.pyplot as plt
from dateutil.parser import parse


# Transponder uplink and downlink frequencies (Hz)
UPLINK_FREQUENCY = 436e6
DOWNLINK_FREQUENCY = 146e6

# Doppler calculation rate.
# The output data from this script will be interpolated when applying it to a test waveform
SAMPLE_RATE = 100

# Pass Information.
# The pass time was computed using gpredict.
# Obsever location
OBSERVER_LAT = -34.43
OBSERVER_LON = 138.72
OBSERVER_ALT = 0.0

# CAS-4A as an example of a linear transponder sat.
TLE_0 = "0 ZHUHAI-1 01"
TLE_1 = "1 42761U 17034D   20087.81104596  .00000284  00000-0  30684-4 0  9995"
TLE_2 = "2 42761  43.0186 114.9591 0012544 242.5178 286.4411 15.09851347153625"

PASS_START = "2020-03-28T05:25Z"
PASS_DURATION = 13*60 # seconds

TIMESTEP = 1/SAMPLE_RATE

observer = ephem.Observer()
observer.lat = str(OBSERVER_LAT)
observer.lon = str(OBSERVER_LON)
observer.elevation = OBSERVER_ALT

satellite = ephem.readtle(TLE_0, TLE_1, TLE_2)

def calculate_velocity(
    tle: list,
    observer_loc: list,
    timestamp
):
    """ Calculate the relative velocity between an observer and a satellite, at a given time """

    observer.date = timestamp.strftime('%Y-%m-%d %H:%M:%S.%f')

    satellite.compute(observer)

    return satellite.range_velocity


def compute_pass():
    """ Compute the relative velocity over a satellite pass """

    _start_time = parse(PASS_START)

    _time = []
    _velocity = []

    for _time_offset in np.arange(0, PASS_DURATION, TIMESTEP):
        _vel = calculate_velocity(
            [TLE_0, TLE_1, TLE_2],
            [OBSERVER_LAT, OBSERVER_LON, OBSERVER_ALT],
            _start_time + datetime.timedelta(0,_time_offset)
        )

        _time.append(_time_offset)
        _velocity.append(_vel)

        if _time_offset%10.0 == 0:
            print("Time: %.1f" % _time_offset)

    return (np.array(_time), np.array(_velocity))


if __name__ == "__main__":

    (_time, _velocity) = compute_pass()

    # Calculate the offset from the uplink centre freq as observer at the satellite.
    _sat_observed_uplink = UPLINK_FREQUENCY * (1 - _velocity/ ephem.c) - UPLINK_FREQUENCY

    # Calculate the 'raw' downlink doppler as observed at the ground station.
    _downlink_doppler_raw = DOWNLINK_FREQUENCY * (1 - _velocity/ ephem.c) - DOWNLINK_FREQUENCY

    # Calculte the observer 
    _downlink_doppler = (-1*_sat_observed_uplink + DOWNLINK_FREQUENCY) * (1 - _velocity/ ephem.c) - DOWNLINK_FREQUENCY

    # Save doppler data to disk
    np.savez('doppler.npz', _time, _downlink_doppler)

    # Plot
    plt.plot(_time, _sat_observed_uplink, label="Observed uplink")
    plt.plot(_time, _downlink_doppler_raw, label="Raw Downlink Doppler")
    plt.plot(_time, _downlink_doppler, label="Observed downlink")
    plt.legend()
    plt.show()

