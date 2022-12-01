#!/usr/bin/env bash
# snr_curves.sh
#
# Plot SNR curves for raw data modes

set -x

PATH=${PATH}:${HOME}/codec2/build_linux/src
CODEC2=${HOME}/codec2

# run over a range of No values
# pull our SNR meas, SNR est, PER
# run all three modes
# run AWGN + MPP channels

No_list='-20 -18 -16 -14 -12 -10'
Ntestframes=10
mode='datac3'

function generate_curve_data {
  mode=$1
  tx_log=$(mktemp)
  ch_log=$(mktemp)
  rx_log=$(mktemp)

  i=1
  rm snrest_${mode}*.txt
  for No in $No_list
  do
    freedv_data_raw_tx --bursts $Ntestframes --testframes $Ntestframes $mode /dev/zero - 2>${tx_log} | \
    ch - - --No $No -f 20 2>>${ch_log} | \
    freedv_data_raw_rx --testframes $mode - /dev/null 2>${rx_log} -v
    SNRoffset=$(cat ${tx_log} | grep "mark:space" | tr -s ' ' | cut -d' ' -f 5)
    echo ${SNRoffset} > snroffset_${mode}.txt
    SNRest=$(cat ${rx_log} | grep '\-BS\-' | tr -s ' ' | cut -d' ' -f17)
    if [ ! -z "$SNRest" ]; then
      echo ${SNRest} > snrest_${mode}_${i}.txt
    fi
    PERmeas=$(cat ${rx_log} | grep 'Coded FER' | cut -d' ' -f3)
    echo ${PERest} >> per_${mode}.txt
    i=$((i+1))
  done

  SNRch=$(cat ${ch_log} | grep SNR3k | tr -s ' ' | cut -d' ' -f3)
  echo ${SNRch} > snrch_${mode}.txt
}

generate_curve_data 'datac0'
generate_curve_data 'datac1'
generate_curve_data 'datac3'
