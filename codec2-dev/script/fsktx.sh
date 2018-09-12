#!/bin/sh
# fsktx.sh
# David Rowe Sep 2018
#
# Transmit 100 baud 2FSK test frames using Rpitx
#
# usage: $ sudo /fsktx.sh FreqHz Bits
#        $ sudo ./fsktx.sh 7177000 6000

CODEC2_BIN=/home/david/codec2-dev/build_linux/src
RPITX_BIN=/home/david/rpitx2

$CODEC2_BIN/fsk_get_test_bits - $2 | $CODEC2_BIN/fsk_mod_ext_vco - $RPITX_BIN/2fsk.f 2 --rpitx 800 100
sudo $RPITX_BIN/freedv $RPITX_BIN/2fsk.f $1 100
