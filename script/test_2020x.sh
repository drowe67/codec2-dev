#!/bin/bash -x
# test_2020x.sh
# David Rowe Feb 2022
#
# Script to support testing experimental 2020A and 2020B modes and 700E control.

CODEC2_PATH=$HOME/codec2
PATH=$PATH:$CODEC2_PATH/build_linux/src:$CODEC2_PATH/build_linux/misc
FADING_DIR=$CODEC2_PATH/build_linux/unittest
No_AWGN=-20
No_Multipath=-25

function run_sim() {
    fullfile=$1
    filename=$(basename -- "$fullfile")
    extension="${filename##*.}"
    filename="${filename%.*}"
    mode=$2
    if [ "$mode" == "700E" ] || [ "$mode" == "700D" ]; then
        rateHz=8000
    else
        rateHz=16000
    fi
    clip=$3
    if [ "$clip" == "clip" ]; then
        clipflag=1
    else
        clipflag=0
    fi       
    channel=$4
    if [ "$channel" == "awgn" ]; then
        channel_opt=""
        No=$No_AWGN
    else
        channel_opt='--'${channel}
        No=$No_Multipath
    fi
    
    freedv_tx ${mode} ${fullfile} - --clip ${clipflag} | \
    ch - - --No $No ${channel_opt} --fading_dir ${FADING_DIR} | \
    freedv_rx ${mode} - - | \
    sox -t .s16 -r ${rateHz} -c 1 - ${filename}_${mode}_${clip}_${channel}.wav
}

run_sim ~/LPCNet/wav/peter.wav 2020 noclip awgn
run_sim ~/codec2/wav/big_dog.wav 700E clip mpp
