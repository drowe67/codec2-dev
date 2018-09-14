#!/bin/sh
# fsktxcots.sh
# David Rowe Sep 2018
#
# Transmit 100 baud 2FSK test frames using COTS radio connected via USB
#
# usage: $ ./fsktxcots.sh Bits
#        $ ./fsktxcots.sh 6000

CODEC2_BIN=/home/david/codec2-dev/build_linux/src

echo 'T 1' | rigctl -m 361 -r /dev/ttyUSB0
$CODEC2_BIN/fsk_get_test_bits - $1 | $CODEC2_BIN/fsk_mod 2 8000 100 1200 100 - - | aplay -f S16_LE
echo 'T 0' | rigctl -m 361 -r /dev/ttyUSB0

