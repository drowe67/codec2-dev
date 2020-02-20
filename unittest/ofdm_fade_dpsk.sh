#!/bin/bash -x
# 
# David Sep 2019
# Tests 2020 OFDM modem fading channel performance in DPSK mode, using a simulated faster (2Hz) high SNR fading channel

PATH=$PATH:../build_linux/src
RAW=$PWD/../raw
results=$(mktemp)

# generate fading file
if [ ! -f ../raw/faster_fading_samples.float ]; then
    echo "Generating fading files ......"
    cmd='cd ../octave; pkg load signal; cohpsk_ch_fading("../raw/faster_fading_samples.float", 8000, 2.0, 8000*60)'
    octave --no-gui -qf --eval "$cmd"
    [ ! $? -eq 0 ] && { echo "octave failed to run correctly .... exiting"; exit 1; }
fi

pwd
# Coded BER should be < 1% for this test
nc=31; ofdm_mod --in /dev/zero --testframes 300 --nc $nc --ldpc 2 --verbose 1 -p 312 --dpsk | cohpsk_ch - - -40 --Fs 8000 -f 10 --ssbfilt 1 --faster --raw_dir $RAW | ofdm_demod --out /dev/null --testframes --nc $nc --verbose 1 --ldpc 2 -p 312 --dpsk 2> $results
cat $results
cber=$(cat $results | sed -n "s/^Coded BER.* \([0-9..]*\) Tbits.*/\1/p")
python -c "import sys; sys.exit(0) if $cber<=0.05 else sys.exit(1)"

