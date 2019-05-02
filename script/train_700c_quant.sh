#!/bin/bash -x
# train_700C_quant.sh
# David Rowe May 2019
#
# Train a Vector Quantiser (VQ) for Codec 2 700C
# This is a two stage VQ with 512 entries (9 bits) per stage

SRC=~/Downloads/all_speech_8k.sw
CODEC2_BUILD=/home/david/codec2/build_linux
K=20

# c2enc can dump "feature vectors" that contain the current VQ input
$CODEC2_BUILD/src/c2enc 700C $SRC /dev/null --mlfeat feat.f32
# extract VQ input as trarining data, thena train two stage VQ
$CODEC2_BUILD/misc/extract -s 1 -e $K -t 41 feat.f32 stage0_in.f32
$CODEC2_BUILD/misc/vqtrain stage0_in.f32 $K 512 vq_stage1.f32 -s 1e-3 -r stage1_in.f32
$CODEC2_BUILD/misc/vqtrain stage1_in.f32 $K 512 vq_stage2.f32 -s 1e-3




