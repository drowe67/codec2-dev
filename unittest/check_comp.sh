#!/usr/bin/env bash
#
# Check Octave and C raw data mode waveforms have about the same
# compression - sanity check for C port of raw data modes.
#
# For manual run outside of ctest:
#   cd codec/build_linux
#  ../unittest/check_comp.sh ${CODEC2} ${PATH}:${CODEC2}/build_linux/src

CODEC2=$1
PATH=$2:$PATH
set -x
octave_log=$(mktemp)
ch_log=$(mktemp)

echo "warning ('off', 'Octave:data-file-in-path');
      ofdm_ldpc_tx('test_datac0.raw','datac0',1,100,'awgn','bursts',10,'txclip'); 
      quit" | DISPLAY="" octave-cli -p ${CODEC2}/octave 1>${octave_log}
oct_rms=$(cat ${octave_log} | tr -s ' ' | grep 'RMS:' | cut -d' ' -f4)
oct_cpapr=$(cat ${octave_log} | grep 'RMS:' | tr -s ' ' | cut -d' ' -f6)
    
freedv_data_raw_tx datac0 /dev/zero - --delay 1000 --testframes 10 --bursts 10 --clip 1 --txbpf 1 | \
ch - /dev/null 2>${ch_log}
ch_rms=$(cat ${ch_log} | grep RMS | tr -s ' ' | cut -d' ' -f5)
ch_cpapr=$(cat ${ch_log} | grep RMS | tr -s ' ' | cut -d' ' -f7)

# Allow 5% difference
python3 -c "import sys; sys.exit(0) if abs((${oct_rms} - ${ch_rms})/${oct_rms}) < 0.05 else sys.exit(1)"
python3 -c "import sys; sys.exit(0) if abs((${oct_cpapr} - ${ch_cpapr})/${oct_cpapr}) < 0.05 else sys.exit(1)"
