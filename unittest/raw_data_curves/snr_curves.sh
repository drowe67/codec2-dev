# snr_curves.sh
#
# Library of bash functions to generate data for SNR curves.
#
# testing a function example:
#   $ bash -c "source ./snr_curves.sh; generate_octave_tx_data datac0 awgn"

set -x

PATH=${PATH}:${HOME}/codec2/build_linux/src
CODEC2=${HOME}/codec2
FADING_DIR=${CODEC2}/build_linux/unittest

snr_list='-5 -4 -3 -2 0 1 2 4'
No_list='-13 -14 -15 -16 -18 -20 -22 -24 -26'
Nbursts_awgn=20
Nbursts_mpp=100

# Octave Tx injects noise and is source of truth for SNR, measure BER/PER v SNR
function generate_octave_tx_data {
  mode=$1
  channel=$2
  
  Nbursts=$Nbursts_awgn
  snr_nudge=0
  if [ "$channel" == "mpp" ]; then
    Nbursts=$Nbursts_mpp
    snr_nudge=4  
  fi
  
  rx_log=$(mktemp)

  i=1
  rm -f snr_oct_${mode}_${channel}*.txt
  rm -f ber_oct_${mode}_${channel}*.txt
  rm -f per_oct_${mode}_${channel}*.txt
  for snr in $snr_list
  do
    snr_adj=$((${snr}+${snr_nudge}))  
    echo "warning ('off', 'Octave:data-file-in-path');
    ofdm_ldpc_tx('test_${mode}.raw','${mode}',1,${snr_adj},'${channel}','bursts',${Nbursts},'crc');
    quit" | DISPLAY="" octave-cli -p ${CODEC2}/octave
    freedv_data_raw_rx --testframes $mode test_${mode}.raw /dev/null 2>${rx_log} -v
    BERmeas=$(cat ${rx_log} | grep 'BER......:' | cut -d' ' -f2)
    PERmeas=$(cat ${rx_log} | grep 'Coded FER' | cut -d' ' -f3)
    
    echo ${snr_adj} >> snr_oct_${mode}_${channel}.txt
    echo ${BERmeas} >> ber_oct_${mode}_${channel}.txt
    echo ${PERmeas} >> per_oct_${mode}_${channel}.txt
    i=$((i+1))
  done
  echo 0 > offset_oct_${mode}_${channel}.txt
}

# ch injects noise and is source of truth for SNR, measure BER/PER v SNR
# Octave Tx
function generate_ch_data {
  mode=$1
  channel=$2

  ch_multipath=''
  Nbursts=$Nbursts_awgn
  snr_nudge=0  
  if [ "$channel" == "mpp" ]; then
    ch_multipath='--mpp'
    Nbursts=$Nbursts_mpp
    snr_nudge=4  
  fi    
    
  octave_log=$(mktemp)
  ch_log=$(mktemp)
  rx_log=$(mktemp)

  i=1
  rm -f snr_ch_${mode}_${channel}*.txt
  rm -f ber_ch_${mode}_${channel}*.txt
  rm -f per_ch_${mode}_${channel}*.txt
  for No in $No_list
  do
    No_adj=$((${No}-${snr_nudge}))  
    echo "warning ('off', 'Octave:data-file-in-path');
    ofdm_ldpc_tx('test_${mode}.raw','${mode}',1,100,'awgn','bursts',${Nbursts},'crc'); 
    quit" | DISPLAY="" octave-cli -p ${CODEC2}/octave 1>${octave_log}
    SNRoffset=$(cat ${octave_log} | grep 'Burst offset:' | cut -d' ' -f5)
    
    ch test_${mode}.raw - --No $No_adj ${ch_multipath} --fading_dir ${FADING_DIR} 2>>${ch_log} | \
    freedv_data_raw_rx --testframes $mode - /dev/null -v 2>${rx_log}
    BERmeas=$(cat ${rx_log} | grep 'BER......:' | cut -d' ' -f2)
    PERmeas=$(cat ${rx_log} | grep 'Coded FER' | cut -d' ' -f3)
    
    echo ${BERmeas} >> ber_ch_${mode}_${channel}.txt
    echo ${PERmeas} >> per_ch_${mode}_${channel}.txt
    i=$((i+1))
  done
  
  echo ${SNRoffset} > offset_ch_${mode}_${channel}.txt
  SNRch=$(cat ${ch_log} | grep SNR3k | tr -s ' ' | cut -d' ' -f3)
  echo ${SNRch} > snr_ch_${mode}_${channel}.txt
}

# ch injects noise and is source of truth for SNR, measure BER/PER v SNR and
# SNR estimates v SNR from rx, C Tx
function generate_snrest_v_snr_data {
  mode=$1
  channel=$2  

  snr_nudge=0
  aNo_list=$No_list

  # nudge SNR test range to get meaningful results for these tests  
  if [ "$mode" == "datac1" ]; then
    snr_nudge=4
  fi
  if [[ "$mode" == "datac4" || "$mode" == "datac13" ]]; then
    snr_nudge=-6
  fi
  
  ch_multipath=''
  Nbursts=$Nbursts_awgn
  if [ "$channel" == "mpp" ]; then
    ch_multipath='--mpp'
    Nbursts=$Nbursts_mpp
    snr_nudge=$((${snr_nudge}+4))
  fi    

  clip=0
  id='ctx'
  if [ "$#" -eq 3 ]; then
    clip=$3
    id='ctxc'
    snr_nudge=$((${snr_nudge}-4))
  fi

  tx_log=$(mktemp)
  ch_log=$(mktemp)
  rx_log=$(mktemp)

  i=1
  rm -f snrest_${id}_${mode}_${channel}*.txt
  rm -f ber_${id}_${mode}_${channel}*.txt
  rm -f per_${id}_${mode}_${channel}*.txt
  for No in $aNo_list
  do
    No_adj=$((${No}-${snr_nudge}))  
    freedv_data_raw_tx --clip ${clip} --delay 1000 --txbpf ${clip} --bursts $Nbursts --testframes $Nbursts $mode /dev/zero - 2>${tx_log} | \
    ch - - --No $No_adj ${ch_multipath} --fading_dir ${FADING_DIR} 2>>${ch_log} | \
    freedv_data_raw_rx --testframes $mode - /dev/null 2>${rx_log} -v
    SNRoffset=$(cat ${tx_log} | grep "mark:space" | tr -s ' ' | cut -d' ' -f 5)
   
    SNRest=$(cat ${rx_log} | grep '\-BS\-' | tr -s ' ' | cut -d' ' -f17)
    if [ ! -z "$SNRest" ]; then
      echo ${SNRest} > snrest_${id}_${mode}_${channel}_${i}.txt
    fi
    BERmeas=$(cat ${rx_log} | grep 'BER......:' | cut -d' ' -f2)
    PERmeas=$(cat ${rx_log} | grep 'Coded FER' | cut -d' ' -f3)
    echo ${BERmeas} >> ber_${id}_${mode}_${channel}.txt
    echo ${PERmeas} >> per_${id}_${mode}_${channel}.txt
    i=$((i+1))
  done

  echo ${SNRoffset} > offset_${id}_${mode}_${channel}.txt
 
  # trap not enough fading file samples (with mpp)
  grep "Fading file finished" ${ch_log}
  if [ $? -eq 0 ]; then
      cat ${ch_log}
      exit 1
  fi
  SNRch=$(cat ${ch_log} | grep SNR3k | tr -s ' ' | cut -d' ' -f3)
  echo ${SNRch} > snr_${id}_${mode}_${channel}.txt
}

# Sanity check to make sure Octave/CML is set up OK
function test_ldpc {
  echo "ldpcut; quit" | DISPLAY="" octave-cli -p ${CODEC2}/octave
  if [ "$?" -ne 0 ]; then
    echo "basic octave test failed, you may need to"
    echo "(a) run ctests to create build_xxx/cml"
    echo "(b) set up ~/.octaverc as per octave/ldpc.m"
    exit 1
  else
      echo "OK"
  fi
}
