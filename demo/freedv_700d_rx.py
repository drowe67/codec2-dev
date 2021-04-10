#!/usr/bin/env python3
'''
  Demo receive program for FreeDV API 700D mode.

  cd ~/codec2/build_linux
  cat ../raw/ve9qrp_10s.raw | ./demo/freedv_700d_tx | ../demo/freedv_700d_rx.py | aplay -f S16_LE

  Credits: Thanks DJ2LS, xssfox, VK5QI
'''

import ctypes
from ctypes import *
import sys
import pathlib

libname = pathlib.Path().absolute() / "src/libcodec2.so"
c_lib = ctypes.CDLL(libname)
c_lib.freedv_open.restype = ctypes.POINTER(ctypes.c_ubyte)
c_lib.freedv_get_n_max_modem_samples.restype = c_int
c_lib.freedv_get_n_max_speech_samples.restype = c_int
c_lib.freedv_nin.restype = c_int

FREEDV_MODE_700D = 7 # from freedv_api.h             
freedv = c_lib.freedv_open(FREEDV_MODE_700D)

n_max_speech_samples = int(c_lib.freedv_get_n_max_speech_samples(freedv))
n_max_modem_samples = c_lib.freedv_get_n_max_modem_samples(freedv)
c_lib.freedv_rawdatarx.argtype = [POINTER(c_ubyte), c_short * n_max_speech_samples, c_short * n_max_modem_samples]

speech_out = (ctypes.c_short * n_max_speech_samples)()
speech_out = bytes(speech_out)

while True:
    nin = c_lib.freedv_nin(freedv)
    demod_in = sys.stdin.buffer.read(nin*2)
    if len(demod_in) == 0: quit()
    nout = c_lib.freedv_rx(freedv, speech_out, demod_in)
    sys.stdout.buffer.write(speech_out[:nout*2])
