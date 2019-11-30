#!/bin/bash -x
# 
# David June 2019
# Tests 700D OFDM modem fading channel performance, using a simulated channel

PATH=$PATH:../build_linux/src
RAW=$PWD/../raw
results=$(mktemp)

# generate fading file
if [ ! -f ../raw/fast_fading_samples.float ]; then
    echo "Generating fading files ......"
    cmd='cd ../octave; pkg load signal; cohpsk_ch_fading("../raw/fast_fading_samples.float", 8000, 1.0, 8000*60)'
    octave --no-gui -qf --eval "$cmd"
    [ ! $? -eq 0 ] && { echo "octave failed to run correctly .... exiting"; exit 1; }
fi

pwd
# BER should be around 4% for this test (it's better for larger interleavers but no one uses interleaving in practice)
ofdm_mod --in /dev/zero --ldpc 1 --testframes 60 --txbpf | cohpsk_ch - - -24 --Fs 8000 -f -10 --fast --raw_dir $RAW | ofdm_demod --out /dev/null --testframes --verbose 2 --ldpc 1 2> $results
cat $results
cber=$(cat $results | sed -n "s/^Coded BER.* \([0-9..]*\) Tbits.*/\1/p")
python -c "import sys; sys.exit(0) if $cber<=0.05 else sys.exit(1)"

