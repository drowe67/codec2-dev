#!/bin/bash -x
# ratek_resampler.sh
# David Rowe Sep 2022
#
# Support for rate K resampler experiments see doc/ratek_resampler

CODEC2_PATH=$HOME/codec2
PATH=$PATH:$CODEC2_PATH/build_linux/src:$CODEC2_PATH/build_linux/misc
K=30
M=512
Kst=0
Ken=29

# train a new VQ
function train() {
  fullfile=$1
  filename=$(basename -- "$fullfile")
  extension="${filename##*.}"
  filename="${filename%.*}"
  
  extract -t $K -s $Kst -e $Ken --removemean --writeall $fullfile ${filename}_nomean.f32
  vqtrain ${filename}_nomean.f32 $K $M  --st $Kst --en $Ken -s 1e-3 vq_stage1.f32 -r res1.f32
  vqtrain res1.f32 $K $M  --st $Kst --en $Ken  -s 1e-3 vq_stage2.f32 -r res2.f32
}

train $1

