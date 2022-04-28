#!/bin/bash
# train_comp.sh
# David Rowe April 2022
#
# Training and testing Vector Quantisers (VQ) for Codec 2 newamp1, using
# audio compression

TRAIN=~/Downloads/train.spc
CODEC2_PATH=$HOME/codec2
PATH=$PATH:$CODEC2_PATH/build_linux/src:$CODEC2_PATH/build_linux/misc
K=30
Kst=2
Ken=16

# compressor with AGC
function compress() {
  fullfile=$1
  filename=$(basename -- "$fullfile")
  extension="${filename##*.}"
  filename="${filename%.*}"
  tmp=$(mktemp)
  # SSB filter -> peak after filter (removes F0)
  ch $fullfile - 2>/dev/null | ch - /dev/null 2>$tmp
  peak=$(cat $tmp | grep peak | tr -s ' ' | cut -f3 -d' ')
  gain=$(python -c "gain=32767.0/${peak}; print(\"%3.1f\" % gain)")
  ch $fullfile ${filename}_comp.s16 --gain $gain --clip 16384 2>$tmp
  papr=$(cat $tmp | grep peak | tr -s ' ' | cut -f7 -d' ')
  printf "PAPR: %5.2f\n" $papr
}

# train a new VQ
function train() {
  fullfile=$TRAIN
  filename=$(basename -- "$fullfile")
  extension="${filename##*.}"
  filename="${filename%.*}"
  
  c2sim $fullfile --rateK --rateKout ${filename}.f32
  echo "ratek=load_f32('../build_linux/${filename}.f32',20); vq_700c_eq; ratek_lim=limit_vec(ratek, 0, 40); save_f32('../build_linux/${filename}_lim.f32', ratek_lim); quit" | \
  octave -p ${CODEC2_PATH}/octave -qf
  vqtrain ${filename}_lim.f32 $K 4096 vq_stage1.f32 -s 1e-3 --st $Kst --en $Ken
}

function listen() {
  fullfile=$1
  filename=$(basename -- "$fullfile")
  extension="${filename##*.}"
  filename="${filename%.*}"
  
  c2sim $fullfile --rateK --rateKout ${filename}.f32
  echo "ratek=load_f32('../build_linux/${filename}.f32',20); vq_700c_eq; ratek_lim=limit_vec(ratek, 0, 40); save_f32('../build_linux/${filename}_lim.f32', ratek_lim); quit" | \
  octave -p ${CODEC2_PATH}/octave -qf
  cat ${filename}_lim.f32 | vq_mbest --st $Kst --en $Ken -k $K -q vq_stage1.f32 > ${filename}_test.f32
  c2sim $fullfile --rateK --rateKin ${filename}_test.f32 -o - | sox -t .s16 -r 8000 -c 1 - ${filename}_sub.wav
  c2sim $fullfile --rateK --newamp1vq -o - | sox -t .s16 -r 8000 -c 1 - ${filename}_newamp1.wav
}

function run() {
    # choose which function to run here
    train
    # these two samples are inside training database
    listen ~/Downloads/fish_8k.sw
    listen ~/Downloads/cap_8k.sw
    # two samples from outside training database
    listen $CODEC2_PATH/raw/big_dog.raw
    listen $CODEC2_PATH/raw/hts2a.raw
    # these two samples are inside training database, but with LPF at 3400 Hz outside of subset
    listen ~/Downloads/fish_8k_lp.sw
    listen ~/Downloads/cap_8k_lp.sw
}

compress $CODEC2_PATH/raw/big_dog.raw
compress $CODEC2_PATH/raw/hts2.raw
compress ~/Downloads/fish.s16
compress ~/Downloads/pencil.s16
