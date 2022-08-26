#!/bin/bash -x
# phase_exp.sh
# David Rowe November 2021
#
# Script to support phase experiments

CODEC2_PATH=$HOME/codec2
PATH=$PATH:$CODEC2_PATH/build_linux/src:$CODEC2_PATH/build_linux/misc

function process() {
  rawfile=$1
  fn=$(basename -- "$rawfile")
  fn="${fn%.*}"
  # 1- original amplitude and phase
  # 2- rateK (coarsely sampled/smoothed K=20) amplitude and orig phase
  # 3- original amplitude and phase0 (derived from amplitude via Hilbert transform) phase
  # 4- rateK amplitude and phase0 (derived from rateK amplitudes)
  # 5- rateK (less coarsely sampled K=30) amplitude and phase0 (derived from rateK=30 amplitudes)
  c2sim ${rawfile} -o - | sox -t .s16 -r 8000 -c 1 - ${fn}_1_orig.wav
  c2sim ${rawfile} --rateK -o - | sox -t .s16 -r 8000 -c 1 - ${fn}_2_ratek.wav
  c2sim ${rawfile} --phase0 --postfilter -o - | sox -t .s16 -r 8000 -c 1 - ${fn}_3_p0.wav
  c2sim ${rawfile} --rateK --phase0 --postfilter -o - | sox -t .s16 -r 8000 -c 1 - ${fn}_4_ratek_p0.wav
  c2sim ${rawfile} --rateK --setK 30 --phase0 --postfilter -o - | sox -t .s16 -r 8000 -c 1 - ${fn}_5_ratek_p0_k30.wav
}

process ${CODEC2_PATH}/raw/hts1a.raw
process ${CODEC2_PATH}/raw/hts2a.raw
process ${CODEC2_PATH}/raw/big_dog.raw
process ${CODEC2_PATH}/raw/g3plx.raw
process ${CODEC2_PATH}/raw/mmt1.raw
process ${CODEC2_PATH}/raw/kristoff.raw
process ${CODEC2_PATH}/raw/cq_ref.raw
process ${CODEC2_PATH}/raw/ve9qrp_10s.raw
