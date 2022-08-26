#!/bin/bash -x
# plphase.sh
# David Rowe Nov 2021
#
# Support for speech PAPR/phase experiments https://github.com/drowe67/codec2/pull/250

TRAIN=~/Downloads/all_speech_8k.sw
CODEC2_PATH=$HOME/codec2
PATH=$PATH:$CODEC2_PATH/build_linux/src:$CODEC2_PATH/build_linux/misc

if [ $# -lt 1 ]; then
    echo
    echo "Generate input files for octave/plphase.m"
    echo
    echo "  usage:"
    echo "    $ cd ~/codec2/build_linux"
    echo "    $ ./plphase.sh rawFile"
    exit
fi

rawFile=$1
b=$(basename "$rawFile" .raw)

# Generate orig and phase0 versions of output wave file to listen to
c2sim ${rawFile} --rateK -o - | sox -t .s16 -r 8000 -c 1 - ${b}_orig.wav
c2sim ${rawFile} --rateK --phase0 --postfilter -o - | sox -t .s16 -r 8000 -c 1 - ${b}_phase0.wav

# Generate input files for octave/plphase.m
c2sim ${rawFile} --rateK --modelout - | ./misc/est_n0 > ${b}_n0.txt
c2sim ${rawFile} --rateK --phase0 --modelout - | ./misc/est_n0 > ${b}_phase0_n0.txt
c2sim ${rawFile} --rateK --phase0 --dump ${b}
