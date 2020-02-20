#!/bin/bash -x
#
# ofdm_phase_est_bw.sh
# David August 2019

# Tests 2020 OFDM modem phase est bandwidth mode option. Locking the
# phase est bandwidth to "high" is useful for High SNR channels with
# fast fading or high phase noise.  In this test we show that with
# high bandwidth phase est mode, the BER is < 5% for the "--faster" (2
# Hz fading) channel model on a fairly high SNR channel.

PATH=$PATH:../build_linux/src
RAW=$PWD/../raw
results=$(mktemp)

# generate fading file
if [ ! -f ../raw/faster_fading_samples.float ]; then
    echo "Generating fading file ......"
    cmd='cd ../octave; pkg load signal; cohpsk_ch_fading("../raw/faster_fading_samples.float", 8000, 2.0, 8000*60)'
    octave --no-gui -qf --eval "$cmd"
    [ ! $? -eq 0 ] && { echo "octave failed to run correctly .... exiting"; exit 1; }
fi

pwd
# BER should be < 5% for this test
nc=37; ofdm_mod --in /dev/zero --testframes 300 --nc $nc --ldpc 2 --verbose 0 | cohpsk_ch - - -40 --Fs 8000 -f 10 --ssbfilt 1 --faster --raw_dir $RAW | ofdm_demod --out /dev/null --testframes --nc $nc --verbose 1 --ldpc 2 --bandwidth 1 2> $results
cat $results
cber=$(cat $results | sed -n "s/^Coded BER.* \([0-9..]*\) Tbits.*/\1/p")
python -c "import sys; sys.exit(0) if $cber<=0.05 else sys.exit(1)"

