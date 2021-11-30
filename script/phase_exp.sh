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
  c2sim ${rawfile} -o - | sox -t .s16 -r 8000 -c 1 - ${fn}_1_orig.wav
  c2sim ${rawfile} --rateK -o - | sox -t .s16 -r 8000 -c 1 - ${fn}_2_ratek.wav
  c2sim ${rawfile} --rateK --phase0 --postfilter -o - | sox -t .s16 -r 8000 -c 1 - ${fn}_3_ratek_p0.wav
  c2sim ${rawfile} --rateK --phase0 --postfilter --postfilter_newamp1 -o - | sox -t .s16 -r 8000 -c 1 - ${fn}_4_ratek_p0_pf.wav
  c2sim ${rawfile} --phase0 --postfilter -o - | sox -t .s16 -r 8000 -c 1 - ${fn}_5_p0.wav
  c2sim ${rawfile} --phase0 --postfilter --dispersion 1 -o - | sox -t .s16 -r 8000 -c 1 - ${fn}_6_p0_disp1.wav
  c2sim ${rawfile} --phase0 --postfilter --dispersion 2 -o - | sox -t .s16 -r 8000 -c 1 - ${fn}_7_p0_disp2.wav
  c2sim ${rawfile} --rateK --phase0 --postfilter --dispersion 1 -o - | sox -t .s16 -r 8000 -c 1 - ${fn}_8_ratek_p0_disp1.wav
  c2sim ${rawfile} --rateK --phase0 --postfilter --dispersion 2 -o - | sox -t .s16 -r 8000 -c 1 - ${fn}_9_ratek_p0_disp2.wav
  c2sim ${rawfile} --dump ${fn}
}

process ${CODEC2_PATH}/raw/hts1a.raw
process ${CODEC2_PATH}/raw/hts2a.raw
process ${CODEC2_PATH}/raw/big_dog.raw
process ${CODEC2_PATH}/raw/g3plx.raw
process ${CODEC2_PATH}/raw/mmt1.raw
process ${CODEC2_PATH}/raw/kristoff.raw
process ${CODEC2_PATH}/raw/cq_ref.raw
