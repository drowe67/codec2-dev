#!/bin/bash
# fskrx.sh
# David Rowe Sep 2018
#
# Receive 100 baud 2FSK test frames using RTL-SDR V3 dongle and nanosdr
#
# usage: $ sudo ./fskrx.sh FreqHz logFileName.json
#        $ sudo ./fskrx.sh 7177000 20180912-1800.json
#
# To stop:
#
#   $ sudo killall nanorx

CODEC2_BIN=/home/david/codec2-dev/build_linux/src
NANORX_BIN=/home/david/nanosdr-0.75/build/nanorx

sudo $NANORX_BIN/nanorx -i rtlsdr -f $1 -m USB --output-rate=48k | $CODEC2_BIN/fsk_demod -f -t 2 48000 100 - /dev/null 2> >(tee -a $2)
