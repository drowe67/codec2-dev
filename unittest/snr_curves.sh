#!/usr/bin/env bash
# snr_curves.sh
#
# The raw data mode modems can estimate channel SNR from the Rx signal. This
# script (and companion Octave snr_curves_plot.m) plots estimated versus actual
# SNR curves for raw data modes

set -x

PATH=${PATH}:${HOME}/codec2/build_linux/src
CODEC2=${HOME}/codec2

# 1. Use Octave Tx built-in channel simulation as source of truth for actual SNR.
# 2. Check the C "ch" tool calibration against Octave source of truth.
# 3. Plot SNR est versus SNR for clean and clipped waveforms
# 4. Plot PER versus SNR
# 5. Further work: plot PER versus peak power to show benefit of clipping, 
#                  extend to MPP channels

snr_list='-4 -2 0 2 4'
No_list='-14 -16 -18 -20 -22'
Nbursts=10

# Using Octave Tx as source of truth for SNR, generate BER/PER v SNR
function generate_octave_tx_data {
  mode=$1

  rx_log=$(mktemp)

  i=1
  rm snr_oct_${mode}*.txt
  rm ber_oct_${mode}*.txt
  rm per_oct_${mode}*.txt
  for snr in $snr_list
  do
    echo "ofdm_ldpc_tx('test_${mode}.raw','${mode}',1,${snr},'awgn','bursts',10,'crc'); quit" | \
    DISPLAY="" octave-cli -p ${CODEC2}/octave
    freedv_data_raw_rx --testframes $mode test_${mode}.raw /dev/null 2>${rx_log} -v
    BERmeas=$(cat ${rx_log} | grep 'BER......:' | cut -d' ' -f2)
    PERmeas=$(cat ${rx_log} | grep 'Coded FER' | cut -d' ' -f3)
    
    echo ${snr} >> snr_oct_${mode}.txt
    echo ${BERmeas} >> ber_oct_${mode}.txt
    echo ${PERmeas} >> per_oct_${mode}.txt
    i=$((i+1))
  done
}

# Using Octave ch as source of truth for SNR, generate BER/PER v SNR
function generate_ch_data {
  mode=$1

  octave_log=$(mktemp)
  ch_log=$(mktemp)
  rx_log=$(mktemp)

  i=1
  rm snr_ch_${mode}*.txt
  rm ber_ch_${mode}*.txt
  rm per_ch_${mode}*.txt
  for No in $No_list
  do
    echo "ofdm_ldpc_tx('test_${mode}.raw','${mode}',1,100,'awgn','bursts',10,'crc'); quit" | \
    DISPLAY="" octave-cli -p ${CODEC2}/octave 1>${octave_log} 
    SNRoffset=$(cat ${octave_log} | grep 'Burst offset:' | cut -d' ' -f5)
    
    ch test_${mode}.raw - --No $No 2>>${ch_log} | \
    freedv_data_raw_rx --testframes $mode - /dev/null -v 2>${rx_log}
    BERmeas=$(cat ${rx_log} | grep 'BER......:' | cut -d' ' -f2)
    PERmeas=$(cat ${rx_log} | grep 'Coded FER' | cut -d' ' -f3)
    
    echo ${BERmeas} >> ber_ch_${mode}.txt
    echo ${PERmeas} >> per_ch_${mode}.txt
    i=$((i+1))
  done
  
  echo ${SNRoffset} > offset_ch_${mode}.txt
  SNRch=$(cat ${ch_log} | grep SNR3k | tr -s ' ' | cut -d' ' -f3)
  echo ${SNRch} > snr_ch_${mode}.txt
}

# Using ch as source of truth for channel SNR, collect SNR estimates from modem
function generate_snrest_v_snr_data {
  mode=$1

  tx_log=$(mktemp)
  ch_log=$(mktemp)
  rx_log=$(mktemp)

  i=1
  rm snrest_${mode}*.txt
  rm per_${mode}*.txt
  for No in $No_list
  do
    freedv_data_raw_tx --bursts $Nbursts --testframes $Nbursts $mode /dev/zero - 2>${tx_log} | \
    ch - - --No $No -f 20 2>>${ch_log} | \
    freedv_data_raw_rx --testframes $mode - /dev/null 2>${rx_log} -v
    SNRoffset=$(cat ${tx_log} | grep "mark:space" | tr -s ' ' | cut -d' ' -f 5)
    echo ${SNRoffset} > snroffset_${mode}.txt
    SNRest=$(cat ${rx_log} | grep '\-BS\-' | tr -s ' ' | cut -d' ' -f17)
    if [ ! -z "$SNRest" ]; then
      echo ${SNRest} > snrest_${mode}_${i}.txt
    fi
    PERmeas=$(cat ${rx_log} | grep 'Coded FER' | cut -d' ' -f3)
    echo ${PERmeas} >> per_${mode}.txt
    i=$((i+1))
  done

  SNRch=$(cat ${ch_log} | grep SNR3k | tr -s ' ' | cut -d' ' -f3)
  echo ${SNRch} > snrch_${mode}.txt
}

#generate_octave_tx_data 'datac0'
generate_ch_data 'datac0'

#generate_snrest_v_snr_data 'datac0'
#generate_snrest_v_snr_data 'datac1'
#generate_snrest_v_snr_data 'datac3'
