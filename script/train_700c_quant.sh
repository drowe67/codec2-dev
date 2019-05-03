#!/bin/bash 
# train_700C_quant.sh
# David Rowe May 2019
#
# Train a Vector Quantiser (VQ) for Codec 2 700C
# This is a two stage VQ with 512 entries (9 bits) per stage

SRC=~/Downloads/all_speech_8k.sw
CODEC2_BUILD=/home/david/codec2/build_linux
K=20
SAMPLES=~/tmp/c2vec_mid

function train() {
  # c2enc can dump "feature vectors" that contain the current VQ input
  $CODEC2_BUILD/src/c2enc 700C $SRC /dev/null --mlfeat feat.f32
  # extract VQ input as trarining data, then train two stage VQ
  $CODEC2_BUILD/misc/extract -s 1 -e $K -t 41 feat.f32 stage0_in.f32
  $CODEC2_BUILD/misc/vqtrain stage0_in.f32 $K 512 vq_stage1.f32 -s 1e-3 -r stage1_in.f32
  $CODEC2_BUILD/misc/vqtrain stage1_in.f32 $K 512 vq_stage2.f32 -s 1e-3
}

function test_a() {
  b=$(basename "$1" .raw)
  $CODEC2_BUILD/src/c2enc 700C $1'.raw' - --var | $CODEC2_BUILD/src/c2dec 700C - - | sox -q -t .s16 -c 1 -r 8000 -b 16  - $SAMPLES/$b'_a.wav'
}

function test_b() {
  b=$(basename "$1" .raw)
  $CODEC2_BUILD/src/c2enc 700C $1'.raw' - --loadcb 1 vq_stage1.f32 --loadcb 2 vq_stage2.f32 --var | \
  $CODEC2_BUILD/src/c2dec 700C - - --loadcb 1 vq_stage1.f32 --loadcb 2 vq_stage2.f32 | sox -q -t .s16 -c 1 -r 8000 -b 16  - $SAMPLES/$b'_b.wav'
}

function feat() {
  b=$(basename "$1" .raw)
  $CODEC2_BUILD/src/c2enc 700C $1'.raw' /dev/null --mlfeat $b.f32
}
function listen() {
  # generate a bunch of test samples for comparsion 
  mkdir -p $SAMPLES
  RAW_FILES="../raw/hts1a ../raw/hts2a ../raw/vk5qi ../raw/cq_ref ../raw/ve9qrp_10s $HOME/Downloads/ma01_01 $HOME/Downloads/c01_01_8k"
  for f in $RAW_FILES
  do
    test_a $f
    test_b $f
  done
}

train
listen
#feat ../raw/hts1a
#feat $HOME/Downloads/ma01_01



