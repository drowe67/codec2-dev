#!/bin/bash -x
# train_comp.sh
# David Rowe April 2022
#
# Training and testing rateK Vector Quantisers (VQ) for Codec 2, using
# audio compression and AGC.
#
# Results May 2022: nm2 > two band comp

TRAIN=~/Downloads/train2.spc
CODEC2_PATH=$HOME/codec2
PATH=$PATH:$CODEC2_PATH/build_linux/src:$CODEC2_PATH/build_linux/misc
K=20
Kst=2
Ken=18

# ---------------------------------------------------------------------------------
# Utility Functions ---------------------------------------------------------------
# ---------------------------------------------------------------------------------

# Just AGC - 200 Hz HPF and normalise peak levels
function agc() {
  fullfile=$1
  filename=$(basename -- "$fullfile")
  extension="${filename##*.}"
  filename="${filename%.*}"
  tmp=$(mktemp)
  cat $fullfile | hpf | ch - /dev/null --ssbfilt 0 2>$tmp
  peak=$(cat $tmp | grep peak | tr -s ' ' | cut -f3 -d' ')
  gain=$(python -c "gain=32767.0/${peak}; print(\"%3.1f\" % gain)")
  cat $fullfile | hpf | ch - ${filename}_agc.s16 --gain $gain --ssbfilt 0 2>$tmp
  papr=$(cat $tmp | grep peak | tr -s ' ' | cut -f7 -d' ')
  printf "PAPR: %5.2f\n" $papr
}

# ---------------------------------------------------------------------------------
# AGC, remove mean based on mean of Kst to Ken (so mean needs to be quantised separately)
# ---------------------------------------------------------------------------------

function train_agc_nm2() {
  fullfile=$TRAIN
  filename=$(basename -- "$fullfile")
  extension="${filename##*.}"
  filename="${filename%.*}"

  agc ${fullfile}
  c2sim ${filename}_agc.s16 --rateK --rateKout ${filename}_agc.f32
  extract ${filename}_agc.f32 ${filename}_agc_nm.f32 -t $K -s $Kst -e $Ken --removemean --writeall --lower 10
  vqtrain ${filename}_agc_nm.f32 $K 2048 vq_stage1.f32 -s 1e-3 --st $Kst --en $Ken -r stage2_in.f32
  vqtrain stage2_in.f32 $K 512 vq_stage2.f32 -s 1e-3 --st $Kst --en $Ken -r stage3_in.f32
  vqtrain stage3_in.f32 $K 512 vq_stage3.f32 -s 1e-3 --st $Kst --en $Ken
}

function listen_agc_nm() {
  fullfile=$1
  filename=$(basename -- "$fullfile")
  extension="${filename##*.}"
  filename="${filename%.*}"

  o=agc_nm3
  mkdir -p $o
  agc ${fullfile}
  c2sim ${filename}_agc.s16 --rateK --rateK_mean_min 0 --rateK_mean_max 60 --rateKnomeanout ${filename}.f32 \
        --K_st  $Kst --K_en $Ken \
        --phase0 --postfilter --dump ${filename} -o - | sox -t .s16 -r 8000 -c 1 - ${o}/${filename}_ratek.wav
  cat ${filename}.f32 | vq_mbest --st $Kst --en $Ken -k $K -q vq_stage1.f32 > ${filename}_vq.f32
  c2sim ${filename}_agc.s16 --rateK --rateK_mean_min 0 --rateK_mean_max 60 --rateKnomeanin ${filename}_vq.f32  \
        --K_st $Kst  --K_en $Ken \
        --phase0 --postfilter --postfilter_newamp1 -o - | sox -t .s16 -r 8000 -c 1 - ${o}/${filename}_agc_vq.wav
  cat ${filename}.f32 | vq_mbest --st $Kst --en $Ken -k $K -q vq_stage1.f32,vq_stage2.f32 --mbest 5 \
      > ${filename}_vq2.f32
  c2sim ${filename}_agc.s16 --rateK --rateK_mean_min 0 --rateK_mean_max 60 --rateKnomeanin ${filename}_vq2.f32  \
        --K_st $Kst --K_en $Ken \
        --phase0 --postfilter --postfilter_newamp1 -o - | sox -t .s16 -r 8000 -c 1 - ${o}/${filename}_agc_vq2.wav
  cat ${filename}.f32 | \
      vq_mbest --st $Kst --en $Ken -k $K -q vq_stage1.f32,vq_stage2.f32,vq_stage3.f32 --mbest 5 \
      > ${filename}_vq3.f32
  c2sim ${filename}_agc.s16 --rateK --rateK_mean_min 0 --rateK_mean_max 60 --rateKnomeanin ${filename}_vq3.f32  \
        --K_st $Kst --K_en $Ken \
        --phase0 --postfilter --postfilter_newamp1 -o - | sox -t .s16 -r 8000 -c 1 - ${o}/${filename}_agc_vq3.wav
  c2sim ${filename}_agc.s16 --rateK --rateK_mean_min 0 --rateK_mean_max 60 --rateKnomeanin ${filename}_vq3.f32  \
        --K_st $Kst --K_en $Ken \
         -o - | sox -t .s16 -r 8000 -c 1 - ${o}/${filename}_agc_vq3_op.wav
  c2sim $fullfile --rateK --newamp1vq \
         --phase0 --postfilter --postfilter_newamp1 -o - | sox -t .s16 -r 8000 -c 1 - ${o}/${filename}_newamp1b.wav
}

# ---------------------------------------------------------------------------------
# AGC, two band freq domain Hilbert compressor
# ---------------------------------------------------------------------------------

function train_compressed_two_band() {
  fullfile=$TRAIN
  filename=$(basename -- "$fullfile")
  extension="${filename##*.}"
  filename="${filename%.*}"

  agc ${fullfile}
  c2sim ${filename}_agc.s16 --rateK --rateKout ${filename}_comp.f32 --comp 75 --comp_gain 10
  extract ${filename}_comp.f32 ${filename}_comp_nm.f32 -t $K -s $Kst -e $Ken --writeall --lower 20 --removemean
  vqtrain ${filename}_comp_nm.f32 $K 2048 vq_stage1.f32 -s 1e-3 --st $Kst --en $Ken -r stage2_in.f32
  vqtrain stage2_in.f32 $K 512 vq_stage2.f32 -s 1e-3 --st $Kst --en $Ken -r stage3_in.f32
  vqtrain stage3_in.f32 $K 512 vq_stage3.f32 -s 1e-3 --st $Kst --en $Ken
}

# listen with agc, two band compressed
function listen_agc_nm_comp() {
  fullfile=$1
  filename=$(basename -- "$fullfile")
  extension="${filename##*.}"
  filename="${filename%.*}"

  o=agc_nm_comp
  mkdir -p $o
  agc ${fullfile}
  c2sim ${filename}_agc.s16 --rateK --comp 75 --comp_gain 10 --rateKnomeanout ${filename}.f32 \
        --K_st $Kst --K_en $Ken \
        --phase0 --postfilter --dump ${filename} -o - | sox -t .s16 -r 8000 -c 1 - ${o}/${filename}_ratek.wav
  cat ${filename}.f32 | vq_mbest --st $Kst --en $Ken -k $K -q vq_stage1.f32 > ${filename}_vq.f32
  c2sim ${filename}_agc.s16 --rateK --comp 75 --comp_gain 10 --rateKnomeanin ${filename}_vq.f32  \
        --K_st $Kst --K_en $Ken \
        --phase0 --postfilter --postfilter_newamp1 -o - | sox -t .s16 -r 8000 -c 1 - ${o}/${filename}_agc_vq.wav
  cat ${filename}.f32 | vq_mbest --st $Kst --en $Ken -k $K -q vq_stage1.f32,vq_stage2.f32 --mbest 5 \
      > ${filename}_vq2.f32
  c2sim ${filename}_agc.s16 --rateK --comp 75 --comp_gain 10 --rateKnomeanin ${filename}_vq2.f32  \
        --K_st $Kst --K_en $Ken \
        --phase0 --postfilter --postfilter_newamp1 -o - | sox -t .s16 -r 8000 -c 1 - ${o}/${filename}_agc_vq2.wav
  cat ${filename}.f32 | \
      vq_mbest --st $Kst --en $Ken -k $K -q vq_stage1.f32,vq_stage2.f32,vq_stage3.f32 --mbest 5 \
      > ${filename}_vq3.f32
  c2sim ${filename}_agc.s16 --rateK --comp 75 --comp_gain 10 --rateKnomeanin ${filename}_vq3.f32  \
        --K_st $Kst --K_en $Ken \
        --phase0 --postfilter --postfilter_newamp1 -o - | sox -t .s16 -r 8000 -c 1 - ${o}/${filename}_agc_vq3.wav
  c2sim ${filename}_agc.s16 --rateK --comp 75 --comp_gain 10 --rateKnomeanin ${filename}_vq3.f32  \
        --K_st $Kst --K_en $Ken \
        -o - | sox -t .s16 -r 8000 -c 1 - ${o}/${filename}_agc_vq3_op.wav
  c2sim $fullfile --rateK --newamp1vq \
        --K_st $Kst --K_en $Ken \
        --phase0 --postfilter --postfilter_newamp1 -o - | sox -t .s16 -r 8000 -c 1 - ${o}/${filename}_newamp1.wav
}

# Look at some stats when running two band freq domain Hilbert compressor
function stats_compressed_two_band() {
  fullfile=$TRAIN
  filename=$(basename -- "$fullfile")
  extension="${filename##*.}"
  filename="${filename%.*}"

  # initial step applies SSB filter
  cat ${fullfile} | hpf | c2sim - --rateK --rateKout ${filename}.f32
  cat ${fullfile} | hpf | c2sim - --rateK --rateKout ${filename}_comp.f32 --comp 75 --comp_gain 10
}

# -----------------------------------------------------------------------
# Exploring the effect of the postfilter with newamp1 and phase0
#  + listened with headphones
#  + PF improves muffled/clicky but adds tonal artefact and AM modulation
# -----------------------------------------------------------------------

function listen_newamp1() {
  fullfile=$1
  filename=$(basename -- "$fullfile")
  extension="${filename##*.}"
  filename="${filename%.*}"

  o=newamp1_pf
  mkdir -p $o
  c2sim $fullfile --rateK -o - | sox -t .s16 -r 8000 -c 1 - ${o}/${filename}_ratek.wav
  c2sim $fullfile --rateK --phase0 --postfilter -o - \
        | sox -t .s16 -r 8000 -c 1 - ${o}/${filename}_ratek_p0.wav
  c2sim $fullfile --rateK --phase0 --postfilter --postfilter_newamp1 -o - \
        | sox -t .s16 -r 8000 -c 1 - ${o}/${filename}_ratek_p0_pf.wav
  c2sim $fullfile --rateK --newamp1vq -o - | sox -t .s16 -r 8000 -c 1 - ${o}/${filename}_newamp1.wav
  c2sim $fullfile --rateK --newamp1vq --phase0 --postfilter -o - \
        | sox -t .s16 -r 8000 -c 1 - ${o}/${filename}_newamp1_p0.wav
  c2sim $fullfile --rateK --newamp1vq --postfilter_newamp1 --phase0 --postfilter -o - \
        | sox -t .s16 -r 8000 -c 1 - ${o}/${filename}_newamp1_p0_pf.wav
}

function listen_test_newamp1() {
    listen_newamp1 $CODEC2_PATH/raw/big_dog.raw
    listen_newamp1 $CODEC2_PATH/raw/hts2a.raw
    listen_newamp1 ~/Downloads/fish.s16
    listen_newamp1 ~/Downloads/pencil.s16
    listen_newamp1 $CODEC2_PATH/raw/kristoff.raw
    listen_newamp1 $CODEC2_PATH/raw/vk5qi.raw
}

function listen_test_agc_nm() {
    listen_agc_nm $CODEC2_PATH/raw/big_dog.raw
    listen_agc_nm $CODEC2_PATH/raw/hts2a.raw
    listen_agc_nm ~/Downloads/fish.s16
    listen_agc_nm ~/Downloads/pencil.s16
    listen_agc_nm $CODEC2_PATH/raw/kristoff.raw
    listen_agc_nm $CODEC2_PATH/raw/vk5qi.raw
    listen_agc_nm ~/Downloads/vk5dgr_testing_8k.wav
}

function listen_test_agc_nm_comp() {
    listen_agc_nm_comp $CODEC2_PATH/raw/big_dog.raw
    listen_agc_nm_comp $CODEC2_PATH/raw/hts2a.raw
    listen_agc_nm_comp ~/Downloads/fish.s16
    listen_agc_nm_comp ~/Downloads/pencil.s16
    listen_agc_nm_comp $CODEC2_PATH/raw/kristoff.raw
    listen_agc_nm_comp $CODEC2_PATH/raw/vk5qi.raw
    listen_agc_nm_comp ~/Downloads/vk5dgr_testing_8k.wav
}


#train_compressed_two_band
#listen_test_agc_nm_comp

train_agc_nm2
listen_test_agc_nm
