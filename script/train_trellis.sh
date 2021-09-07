#!/bin/bash -x
# train_trellis.sh
# David Rowe August 2021
#
# Training Trellis Vector Quantiser for Codec 2 newamp1, supports octave/trellis.m

TRAIN=~/Downloads/all_speech_8k.sw
CODEC2_PATH=$HOME/codec2
PATH=$PATH:$CODEC2_PATH/build_linux/src:$CODEC2_PATH/build_linux/misc
K=20
Kst=2
Ken=16

# train a new VQ and generate quantised training material
function train() {
  fullfile=$TRAIN
  filename=$(basename -- "$fullfile")
  extension="${filename##*.}"
  filename="${filename%.*}"
  
  c2sim $fullfile --rateK --rateKout ${filename}.f32
  echo "ratek=load_f32('../build_linux/${filename}.f32',20); vq_700c_eq; ratek_lim=limit_vec(ratek, 0, 40); save_f32('../build_linux/${filename}_lim.f32', ratek_lim); quit" | \
  octave -p ${CODEC2_PATH}/octave -qf
  #vqtrain ${filename}_lim.f32 $K 4096 vq_stage1.f32 -s 1e-3 --st $Kst --en $Ken
  cat ${filename}_lim.f32 | vq_mbest --st $Kst --en $Ken -k $K -q vq_stage1.f32 > ${filename}_test.f32
}

# choose which function to run here
train
