#!/bin/bash
# train_trellis.sh
# David Rowe August 2021
#
# Script to support:
# 1. Training Trellis Vector Quantiser for Codec 2 newamp1, supports octave/trellis.m
# 2. VQ sorting/optimisation experiments, octave/vq_compare.m

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
  vqtrain ${filename}_lim.f32 $K 4096 vq_stage1.f32 -s 1e-3 --st $Kst --en $Ken

  # VQ the training file
  cat ${filename}_lim.f32 | vq_mbest --st $Kst --en $Ken -k $K -q vq_stage1.f32 > ${filename}_test.f32
}

function listen() {
  fullfile=$1
  filename=$(basename -- "$fullfile")
  extension="${filename##*.}"
  filename="${filename%.*}"

  fullfile_out=$2
  vq_fn=$3
  EbNodB=$4
  
  sox $fullfile -t raw - | c2sim - --rateK --rateKout ${filename}.f32

  echo "ratek=load_f32('../build_linux/${filename}.f32',20); vq_700c_eq; ratek_lim=limit_vec(ratek, 0, 40); save_f32('../build_linux/${filename}_lim.f32', ratek_lim); quit" | \
  octave -p ${CODEC2_PATH}/octave -qf

  echo "pkg load statistics; vq_compare(action='vq_file', '${vq_fn}', EbNodB=${EbNodB}, '${filename}_lim.f32', '${filename}_test.f32'); quit" \ |
  octave -p ${CODEC2_PATH}/octave -qf
      
  sox $fullfile -t raw - | c2sim - --rateK --rateKin ${filename}_test.f32 -o - | sox -t .s16 -r 8000 -c 1 - ${fullfile_out}
}

function print_help {
    echo
    echo "Trellis/VQ optimisation support script"
    echo
    echo "  usage ./train_trellis.sh [-d] [-t] [-v in.wav out.wav vq.f32 EbNodB]"
    echo
    echo "    -d        debug mode; trace script execution"
    echo "    -t        train VQ and generate a fully quantised version of training vectors"
    echo "    -v        synthesis an output file out.wav from in.raw, using the VQ vq.f32"
    echo
    exit
}

# command line arguments to select function

if [ $# -lt 1 ]; then
    print_help
fi

do_train=0
do_vq=0
POSITIONAL=()
while [[ $# -gt 0 ]]
do
key="$1"
case $key in
    -d)
        set -x	
        shift
    ;;
    -t)
        do_train=1
	shift
    ;;
    -v)
        do_vq=1
	in_wav="$2"
	out_wav="$3"
	vq_fn="$4"
	EbNodB="$5"
        shift
	shift
	shift
	shift
	shift
    ;;
    -h)
        print_help	
    ;;
    *)
    POSITIONAL+=("$1") # save it in an array for later
    shift
    ;;
esac
done
set -- "${POSITIONAL[@]}" # restore positional parameters

if [ $do_train -eq 1 ]; then
    train
fi

if [ $do_vq -eq 1 ]; then
  listen ${in_wav} ${out_wav} ${vq_fn} ${EbNodB}
fi
