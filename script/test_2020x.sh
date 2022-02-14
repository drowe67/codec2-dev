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
serial=0

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
        clip_html="yes"
    else
        clipflag=0
        clip_html="no"
    fi       
    channel=$4
    No=-100
    if [ "$channel" == "awgn" ]; then
        channel_opt=""
        No=$No_AWGN
    fi
    if [ "$channel" == "mpp" ] || [ "$channel" == "mpd" ]; then
        channel_opt='--'${channel}
        No=$No_Multipath
    fi

    indopt=$5
    indopt_flag=""
    indopt_html="no"
    indopt_str=""
    if [ "$indopt" == "indopt" ]; then
        indopt_flag="--indopt 1"
        indopt_str="_indopt"
        indopt_html="yes"
    fi
    if [ "$indopt" == "no_indopt" ]; then
        indopt_flag="--indopt 0"
        indopt_str="_no_indopt"
    fi
    
    fn=${filename}_${mode}_${clip}_${channel}${indopt_str}.wav
    freedv_tx ${mode} ${fullfile} - --clip ${clipflag} ${indopt_flag} | \
    ch - - --No $No ${channel_opt} --fading_dir ${FADING_DIR} | \
    freedv_rx ${mode} - - | \
    sox -t .s16 -r ${rateHz} -c 1 - ${fn} trim 0 6

    echo "<tr>"
    echo "<td><a href=\"${fn}\">${serial}</a></td><td>${mode}</td><td>${clip_html}</td><td>${indopt_html}</td><td>${channel}</td>"
    echo "</tr>"
    serial=$((serial+1))
}

# convert speech input file to format we need
SPEECH_IN_16k_WAV=~/Downloads/speech_orig_16k.wav 
SPEECH_IN_16k_RAW=speech_orig_16k.raw
SPEECH_IN_8k_RAW=speech_orig_8k.raw
sox $SPEECH_IN_16k_WAV -t .s16 $SPEECH_IN_16k_RAW
sox $SPEECH_IN_16k_WAV -t .s16 -r 8000 $SPEECH_IN_8k_RAW

echo "<html><table>"
echo "<tr><th>Serial</th><th>Mode</th><th>Clip</th><th>index_opt</th><th>Channel</th></tr>"

# run simulations

run_sim $SPEECH_IN_16k_RAW 2020 noclip clean
run_sim $SPEECH_IN_8k_RAW 700E clip clean

run_sim $SPEECH_IN_16k_RAW 2020 noclip awgn
run_sim $SPEECH_IN_16k_RAW 2020 noclip mpp
run_sim $SPEECH_IN_16k_RAW 2020 noclip mpd
run_sim $SPEECH_IN_16k_RAW 2020 clip awgn
run_sim $SPEECH_IN_16k_RAW 2020 clip mpp
run_sim $SPEECH_IN_16k_RAW 2020 clip mpd

run_sim $SPEECH_IN_16k_RAW 2020A clip awgn indopt
run_sim $SPEECH_IN_16k_RAW 2020A clip mpp  indopt
run_sim $SPEECH_IN_16k_RAW 2020A clip mpp  no_indopt
run_sim $SPEECH_IN_16k_RAW 2020A clip mpd  indopt
run_sim $SPEECH_IN_16k_RAW 2020A clip mpd  no_indopt

run_sim $SPEECH_IN_16k_RAW 2020B clip awgn indopt
run_sim $SPEECH_IN_16k_RAW 2020B clip mpp  indopt
run_sim $SPEECH_IN_16k_RAW 2020B clip mpp  no_indopt
run_sim $SPEECH_IN_16k_RAW 2020B clip mpd  indopt
run_sim $SPEECH_IN_16k_RAW 2020B clip mpd  no_indopt

run_sim $SPEECH_IN_8k_RAW 700E clip awgn
run_sim $SPEECH_IN_8k_RAW 700E clip mpp
run_sim $SPEECH_IN_8k_RAW 700E clip mpd


exit
