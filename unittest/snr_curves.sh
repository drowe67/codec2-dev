#!/usr/bin/env bash
# snr_curves.sh
#
# The raw data mode modems can estimate channel SNR from the Rx signal. This
# script (and companion Octave snr_curves_plot.m) plots estimated versus actual
# SNR curves for raw data modes

#set -x

PATH=${PATH}:${HOME}/codec2/build_linux/src
CODEC2=${HOME}/codec2

# 1. Use Octave Tx built-in channel simulation as source of truth for actual SNR.
# 2. Check the C "ch" tool calibration against Octave source of truth, using PER/BER.
# 3. Plot SNR est versus SNR for clean and clipped waveforms
# 4. Plot PER versus SNR
# 5. Further work: plot PER versus peak power to show benefit of clipping, 
#                  extend to MPP channels

snr_list='-5 -4 -3 -2 0 1 2 4'
No_list='-13 -14 -15 -16 -18 -20 -22'
Nbursts=20

# Using Octave Tx as source of truth for SNR, generate BER/PER v SNR, uses
# Octave Tx
function generate_octave_tx_data {
  mode=$1

  rx_log=$(mktemp)

  i=1
  rm snr_oct_${mode}*.txt
  rm ber_oct_${mode}*.txt
  rm per_oct_${mode}*.txt
  for snr in $snr_list
  do
    echo "ofdm_ldpc_tx('test_${mode}.raw','${mode}',1,${snr},'awgn','bursts',${Nbursts},'crc'); quit" | \
    DISPLAY="" octave-cli -p ${CODEC2}/octave
    freedv_data_raw_rx --testframes $mode test_${mode}.raw /dev/null 2>${rx_log} -v
    BERmeas=$(cat ${rx_log} | grep 'BER......:' | cut -d' ' -f2)
    PERmeas=$(cat ${rx_log} | grep 'Coded FER' | cut -d' ' -f3)
    
    echo ${snr} >> snr_oct_${mode}.txt
    echo ${BERmeas} >> ber_oct_${mode}.txt
    echo ${PERmeas} >> per_oct_${mode}.txt
    i=$((i+1))
  done
  echo 0 > offset_oct_${mode}.txt
}

# Using Octave ch as source of truth for SNR, generate BER/PER v SNR, uses
# Octave Tx
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
    echo "ofdm_ldpc_tx('test_${mode}.raw','${mode}',1,100,'awgn','bursts',${Nbursts},'crc'); quit" | \
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
# Uses C Tx
function generate_snrest_v_snr_data {
  mode=$1
  clip=0
  id='ctx'
  if [ "$#" -eq 2 ]; then
    clip=$2
    id='ctxc'
  fi

  tx_log=$(mktemp)
  ch_log=$(mktemp)
  rx_log=$(mktemp)

  i=1
  rm snrest_${id}_${mode}*.txt
  rm ber_${id}_${mode}*.txt
  rm per_${id}_${mode}*.txt
  for No in $No_list
  do
    freedv_data_raw_tx --clip ${clip} --txbpf ${clip} --bursts $Nbursts --testframes $Nbursts $mode /dev/zero - 2>${tx_log} | \
    ch - - --No $No -f 20 2>>${ch_log} | \
    freedv_data_raw_rx --testframes $mode - /dev/null 2>${rx_log} -v
    SNRoffset=$(cat ${tx_log} | grep "mark:space" | tr -s ' ' | cut -d' ' -f 5)
   
    SNRest=$(cat ${rx_log} | grep '\-BS\-' | tr -s ' ' | cut -d' ' -f17)
    if [ ! -z "$SNRest" ]; then
      echo ${SNRest} > snrest_${id}_${mode}_${i}.txt
    fi
    BERmeas=$(cat ${rx_log} | grep 'BER......:' | cut -d' ' -f2)
    PERmeas=$(cat ${rx_log} | grep 'Coded FER' | cut -d' ' -f3)
    echo ${BERmeas} >> ber_${id}_${mode}.txt
    echo ${PERmeas} >> per_${id}_${mode}.txt
    i=$((i+1))
  done

  echo ${SNRoffset} > offset_${id}_${mode}.txt
  SNRch=$(cat ${ch_log} | grep SNR3k | tr -s ' ' | cut -d' ' -f3)
  echo ${SNRch} > snr_${id}_${mode}.txt
}


# Sanity check to make sure Octave/CML is set up OK
echo "ldpcut; quit" | DISPLAY="" octave-cli -p ${CODEC2}/octave
if [ "$?" -ne 0 ]; then
    echo "basic octave test failed, you may need to"
    echo "(a) run ctests to create build_xxx/cml"
    echo "(b) set up ~/.octaverc as per octave/ldpc.m"
fi

# These results can be rendered with snr_curves_plot.m

# Compare Octave Tx and ch as SNR source of truth, PER/BER curves
# should be on top of each other
generate_octave_tx_data 'datac0'
generate_ch_data 'datac0'
generate_octave_tx_data 'datac1'
generate_ch_data 'datac1'
generate_octave_tx_data 'datac3'
generate_ch_data 'datac3'

# (a) PER/BER for C TX with & without compression 
# (b) Measure SNR estimates v actual SNR
generate_snrest_v_snr_data 'datac0'
generate_snrest_v_snr_data 'datac1'
generate_snrest_v_snr_data 'datac3'
generate_snrest_v_snr_data 'datac0' 1
generate_snrest_v_snr_data 'datac3' 1

