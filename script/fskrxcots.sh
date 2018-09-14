#!/bin/bash
# fskrxcots.sh
# David Rowe Sep 2018
#
# Receive 100 baud 2FSK test frames using COTS HF radio connected via sound card
#
# usage: $ fskrxcots.sh logFileName.json
#        $ fskrxcots.sh 20180912-1800.json
#
# To stop:
#
#   $ sudo killall arecord

CODEC2_BIN=/home/david/codec2-dev/build_linux/src

arecord -D hw:1,0 -f S16_LE -r 48000 - | $CODEC2_BIN/fsk_demod -f -t 2 48000 100 - /dev/null 2> >(tee -a $1)

