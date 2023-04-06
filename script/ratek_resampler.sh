#!/bin/bash -x
# ratek_resampler.sh
# David Rowe Sep 2022
#
# Support for rate K resampler experiments see doc/ratek_resampler

CODEC2_PATH=$HOME/codec2
PATH=$PATH:$CODEC2_PATH/build_linux/src:$CODEC2_PATH/build_linux/misc

# bunch of options we can set via variables
K="${K:-30}"
M="${M:-4096}"
Kst="${Kst:-0}"  # first index 
Ksp="${Ksp:-9}"  # last element of first vector in split
Ken="${Ken:-29}" # last index max K-1
out_dir="${out_dir:-ratek_out}"
extract_options="${extract_options:-}"
options="${options:-}"
mbest="${mbest:-no}"
removemean="${removemean:---removemean}"
lower=${lower:-10}
meanl2=${meanl2:-}
dr=${dr:-100}
drlate=${drlate:-}
stage2="${stage2:-yes}"
stage3="${stage3:-no}"
Nb=20

function batch_process {
  fullfile=$1
  filename=$(basename -- "$fullfile")
  filename="${filename%.*}"
  batch_opt=$2
  outname=$3
  c2sim_opt=$4
  tmp=$(mktemp)
  
  echo "ratek3_batch;" \
       "ratek3_batch_tool(\"${filename}\", "\
                          "'A_out',\"${filename}_a.f32\"," \
                          "'H_out',\"${filename}_h.f32\"," \
                          "${batch_opt},'logfn',\"${tmp}\"); quit;" \
  | octave -p ${CODEC2_PATH}/octave -qf
  c2sim $fullfile --hpf --phase0 --postfilter --amread ${filename}_a.f32 --hmread ${filename}_h.f32 -o - \
  ${c2sim_opt} | sox -t .s16 -r 8000 -c 1 - ${out_dir}/${filename}_${outname}.wav
  printf "%-10s %-20s %4.2f\n" ${filename} ${outname} $(cat ${tmp}) >> ${out_dir}/zlog.txt
}

# 230331: Building up fully quantised candidates
function vq_test_230331() {
  fullfile=$1
  filename=$(basename -- "$fullfile")
  filename="${filename%.*}"
  extension="${filename##*.}"
  mkdir -p $out_dir
  
  c2sim $fullfile --hpf --modelout ${filename}_model.bin

  # orig amp and phase
  c2sim $fullfile --hpf --modelout ${filename}_model.bin -o - | \
  sox -t .s16 -r 8000 -c 1 - ${out_dir}/${filename}_1_out.wav
 
  # Amps Nb filtered, phase0, rate K=20 resampling, phase postfilter,
  # rate L amp postfilter, normalise energy
  batch_process $fullfile "'K',20,'amp_pf','phase_pf','norm_en'" "2_k20"  

  # 3 x 9 VQ, decimate by 3 (target 1200 bits/s when side info added)
  batch_process $fullfile "'K',20,'amp_pf','phase_pf','norm_en','mic_eq',2,'dec',3, \
  'vq1','../build_linux/train_three_vq1.f32', \
  'vq2','../build_linux/train_three_vq2.f32', \ 
  'vq3','../build_linux/train_three_vq3.f32'"  "3_vq3_27_dec3"

  # 3 x 9 VQ, decimate by 3, 3 bit mean quant (target 1200 bits/s when side info added)
  batch_process $fullfile "'K',20,'amp_pf','phase_pf','norm_en','mic_eq',2,'dec',3,'quant_mean3', \
  'vq1','../build_linux/train_three_vq1.f32', \
  'vq2','../build_linux/train_three_vq2.f32', \ 
  'vq3','../build_linux/train_three_vq3.f32'"  "4_vq3_27_dec3_q3"

  # 1 x 12 VQ, decimate by 3 (target 700 bits/s when side info added)
  batch_process $fullfile "'K',20,'amp_pf','phase_pf','norm_en','mic_eq',2,'dec',3, \
  'vq1','../build_linux/train_k20_vq1.f32'" "5_vq1_12_dec3"

  # 1 x 12 VQ, decimate by 3, 3 bit mean quant (target 700 bits/s when side info added)
  batch_process $fullfile "'K',20,'amp_pf','phase_pf','norm_en','mic_eq',2,'dec',3,'quant_mean3', \
  'vq1','../build_linux/train_k20_vq1.f32'" "6_vq1_12_dec3_q3"

  cat $fullfile | hpf | c2enc 3200 - - | c2dec 3200 - - | sox -t .s16 -r 8000 -c 1 - ${out_dir}/${filename}_7_3200.wav
  
}

# 230323: compressor, mean limiting and quantisation
function comp_test_230323() {
  fullfile=$1
  filename=$(basename -- "$fullfile")
  filename="${filename%.*}"
  mkdir -p $out_dir

  # orig amp and phase
  c2sim $fullfile --hpf --modelout ${filename}_model.bin -o - | \
  sox -t .s16 -r 8000 -c 1 - ${out_dir}/${filename}_1_out.wav
 
  # Amps Nb filtered, phase0, rate K=20 resampling, phase postfilter,
  # rate L amp postfilter, normalise energy
  batch_process $fullfile "'K',20,'amp_pf','phase_pf','norm_en'" "2_k20"  

  # with norm_energy and compression applied to rate K vectors
  batch_process $fullfile "'K',20,'amp_pf','phase_pf','norm_en', 'compress_en'" "3_comp_K"

  # with norm_energy and hilbert compression applied at synthesis
  batch_process $fullfile "'K',20,'amp_pf','phase_pf','norm_en'" "4_comp_hc" "--gainoutlin 3.0"

  # (4) plus mean limiting
  batch_process $fullfile "'K',20,'amp_pf','phase_pf','norm_en','limit_mean'" "5_hc_lim" "--gainoutlin 3.0"

  # (4) with 3 bit mean quantisation
  batch_process $fullfile "'K',20,'amp_pf','phase_pf','norm_en','quant_mean3'" "6_hc_q3" "--gainoutlin 3.0"
}

# 230226: debugging clicks
function vq_test_230226() {
  fullfile=$1
  filename=$(basename -- "$fullfile")
  filename="${filename%.*}"
  mkdir -p $out_dir

  # orig amp and phase
  c2sim $fullfile --hpf --modelout ${filename}_model.bin -o - | \
  sox -t .s16 -r 8000 -c 1 - ${out_dir}/${filename}_1_out.wav
 
  # Amps Nb filtered, phase0, rate K=20 resampling, phase postfilter,
  # rate L amp postfilter
  batch_process $fullfile "'K',20,'amp_pf','phase_pf'" "2_k20"  

  # dec 3, with and without norm_energy
  batch_process $fullfile "'K',20,'amp_pf','phase_pf','dec',3" "3_dec3"
  batch_process $fullfile "'K',20,'amp_pf','phase_pf','dec',3,'norm_en'" "4_dec3_norm"

  # 11 stage 12 bit VQ, dec 3, EQ2, with and without norm_energy
  batch_process $fullfile "'K',20,'amp_pf','phase_pf','dec',3,'mic_eq',2, \
  'vq1','../build_linux/train_k20_vq1.f32'" "5_vq1_dec3"
  batch_process $fullfile "'K',20,'amp_pf','phase_pf','dec',3,'mic_eq',2,'norm_en', \
  'vq1','../build_linux/train_k20_vq1.f32'" "6_vq1_dec3_norm"
}

# 230213: Mic EQ versions 1 & 2
function vq_test_230217() {
  fullfile=$1
  filename=$(basename -- "$fullfile")
  filename="${filename%.*}"
  mkdir -p $out_dir
 
  c2sim $fullfile --hpf --modelout ${filename}_model.bin

  # (1) Amps Nb filtered, phase0, rate K=20 resampling, phase postfilter,
  # rate L amp postfilter, pre-emp
  batch_process $fullfile "'K',20,'amp_pf','phase_pf','pre'" "1_k20"

  # with mic EQ 1
  batch_process $fullfile "'K',20,'amp_pf','phase_pf', \
  'vq1','../build_linux/train_k20_vq1.f32', \
  'vq_en',0,'mic_eq',1,'plot_mic_eq'" "2_k20_eq1"

  # with mic EQ 2
  batch_process $fullfile "'K',20,'amp_pf','phase_pf', \
  'vq1','../build_linux/train_k20_vq1.f32', \
  'vq_en',0,'mic_eq',2,'plot_mic_eq'" "3_k20_eq2"

  # 1 x 12 VQ vanilla
  batch_process $fullfile "'K',20,'amp_pf','phase_pf', \
  'vq1','../build_linux/train_k20_vq1.f32'" "4_k20_vq1"

  # 2 x 12 VQ vanilla
  batch_process $fullfile "'K',20,'amp_pf','phase_pf', \
  'vq1','../build_linux/train_k20_vq1.f32', 
  'vq2','../build_linux/train_k20_vq2.f32'" "5_k20_vq2"

  # 1 x 12 VQ with mic EQ 1
  batch_process $fullfile "'K',20,'amp_pf','phase_pf', \
  'vq1','../build_linux/train_k20_vq1.f32', \
  'mic_eq',1" "6_k20_vq1_eq1"

  # 2 x 12 VQ with mic EQ 1
  batch_process $fullfile "'K',20,'amp_pf','phase_pf', \
  'vq1','../build_linux/train_k20_vq1.f32', \
  'vq2','../build_linux/train_k20_vq2.f32', \
  'mic_eq',1" "7_k20_vq2_eq1"
  
  # 1 x 12 VQ with mic EQ 2
  batch_process $fullfile "'K',20,'amp_pf','phase_pf', \
  'vq1','../build_linux/train_k20_vq1.f32', \
  'mic_eq',2" "8_k20_vq1_eq2"

  # 2 x 12 VQ with mic EQ 2
  batch_process $fullfile "'K',20,'amp_pf','phase_pf', \
  'vq1','../build_linux/train_k20_vq1.f32', \
  'vq2','../build_linux/train_k20_vq2.f32', \
  'mic_eq',2" "9_k20_vq2_eq2"

   cat $fullfile | hpf | c2enc 3200 - - | c2dec 3200 - - | sox -t .s16 -r 8000 -c 1 - ${out_dir}/${filename}_10_3200.wav
}

# 230204: Process sample different VQ designs 1x12, 2x12, 2x9, 3x9,dec2 & 3
function vq_test_230204() {
  fullfile=$1
  filename=$(basename -- "$fullfile")
  filename="${filename%.*}"
  extension="${filename##*.}"
  mkdir -p $out_dir
  
  c2sim $fullfile --hpf --modelout ${filename}_model.bin

  # (1) Amps Nb filtered, phase0, rate K=20 resampling, phase postfilter,
  # rate L amp postfilter, pre-emp, EQ2
  batch_process $fullfile "'K',20,'amp_pf','phase_pf', \
  'vq1','../build_linux/train_k20_vq1.f32', \
  'vq_en',0,'mic_eq',2,'plot_mic_eq'" "1_k20"
 
  # 1 x 12 VQ
  batch_process $fullfile "'K',20,'amp_pf','phase_pf','mic_eq',2, \
  'vq1','../build_linux/train_k20_vq1.f32'" "2_k20_vq1_12"

  # 2 x 12 VQ
  batch_process $fullfile "'K',20,'amp_pf','phase_pf','mic_eq',2, \
  'vq1','../build_linux/train_k20_vq1.f32', \
  'vq2','../build_linux/train_k20_vq2.f32'" "3_k20_vq2_24"

  # 2 x 9 VQ
  batch_process $fullfile "'K',20,'amp_pf','phase_pf','mic_eq',2, \
  'vq1','../build_linux/train_three_vq1.f32',
  'vq2','../build_linux/train_three_vq2.f32'" "4_three_vq2_18"

  # 3 x 9 VQ
  batch_process $fullfile "'K',20,'amp_pf','phase_pf','mic_eq',2, \
  'vq1','../build_linux/train_three_vq1.f32', \
  'vq2','../build_linux/train_three_vq2.f32', \ 
  'vq3','../build_linux/train_three_vq3.f32'"  "5_three_vq3_27"

  # 3 x 9 VQ, decimate by 2 (target 2400 bits/s when side info added)
  batch_process $fullfile "'K',20,'amp_pf','phase_pf','mic_eq',2,'dec', 2, \
  'vq1','../build_linux/train_three_vq1.f32', \
  'vq2','../build_linux/train_three_vq2.f32', \ 
  'vq3','../build_linux/train_three_vq3.f32'"  "6_three_vq3_27_dec2"

  # 3 x 9 VQ, decimate by 3 (target 1200 bits/s when side info added)
  batch_process $fullfile "'K',20,'amp_pf','phase_pf','mic_eq',2,'dec', 3, \
  'vq1','../build_linux/train_three_vq1.f32', \
  'vq2','../build_linux/train_three_vq2.f32', \ 
  'vq3','../build_linux/train_three_vq3.f32'"  "7_three_vq3_27_dec3"

  # 1 x 12 VQ, decimate by 3 (target 700 bits/s)
  batch_process $fullfile "'K',20,'amp_pf','phase_pf','mic_eq',2,'dec', 3, \
  'vq1','../build_linux/train_k20_vq1.f32'" "8_k20_vq1_12_dec3"

  # 2 x 12 VQ, dec by 4 (target 700 bits/s)
  batch_process $fullfile "'K',20,'amp_pf','phase_pf','mic_eq',2,'dec', 4, \
  'vq1','../build_linux/train_k20_vq1.f32', \
  'vq2','../build_linux/train_k20_vq2.f32'" "9_k20_vq2_24_dec4"

   cat $fullfile | hpf | c2enc 3200 - - | c2dec 3200 - - | sox -t .s16 -r 8000 -c 1 - ${out_dir}/${filename}_10_3200.wav

}

# 230204: Dynamic range reduction test, see how it sounds with and without VQ
function dr_vq_test_230204() {
  fullfile=$1
  filename=$(basename -- "$fullfile")
  filename="${filename%.*}"
  mkdir -p $out_dir
  
  c2sim $fullfile --hpf --modelout ${filename}_model.bin
  
  # (1) Amps Nb filtered, phase0, rate K=20 resampling, phase postfilter,
  # rate L amp postfilter, pre-emp
  
  batch_process $fullfile "'K',20,'amp_pf','phase_pf','pre'" "1_k20"

  # (2) with 30dB dynamic range limit
  batch_process $fullfile "'K',20,'amp_pf','phase_pf','pre','DR',30" "2_dr"

  # non pre-emp 1 stage VQ
  batch_process $fullfile "'K',20,'amp_pf','phase_pf', \
  'vq1','../build_linux/train_k20_vq1.f32'" "3_k20_vq1"

  # non pre-emp 2 stage VQ
  batch_process $fullfile "'K',20,'amp_pf','phase_pf', \
  'vq1','../build_linux/train_k20_vq1.f32', \
  'vq2','../build_linux/train_k20_vq2.f32'" "4_k20_vq2"

  # pre-emp 1 stage VQ
  batch_process $fullfile "'K',20,'amp_pf','phase_pf','pre', \
  'vq1','../build_linux/train_pre_vq1.f32'" "5_pre_vq1"

  # pre-emp 2 stage VQ
  batch_process $fullfile "'K',20,'amp_pf','phase_pf','pre', \
  'vq1','../build_linux/train_pre_vq1.f32', \
  'vq2','../build_linux/train_pre_vq2.f32'" "6_pre_vq2"

  # pre-emp, dr limit 1 stage VQ
  batch_process $fullfile "'K',20,'amp_pf','phase_pf','pre','DR',30, \
  'vq1','../build_linux/train_comp_vq1.f32'" "7_dr_vq1"
 
  # pre-emp, dr limit 2 stage VQ
  batch_process $fullfile "'K',20,'amp_pf','phase_pf','pre','DR',30, \
  'vq1','../build_linux/train_comp_vq1.f32', \
  'vq2','../build_linux/train_comp_vq2.f32'" "8_dr_vq2"

}

# 230202: Process samples with simulated quant noise, 1 & 2 stage VQ, weighted search, 
# postfilter at rate L after VQ
function vq_test_230202() {
  fullfile=$1
  filename=$(basename -- "$fullfile")
  filename="${filename%.*}"
  mkdir -p $out_dir
  
  c2sim $fullfile --hpf --modelout ${filename}_model.bin
  
  # 1/ Amps Nb filtered, phase0, phase postfilter, rate Lhigh amp postfilter
  batch_process $fullfile "'amp_pf','phase_pf'" "1_lhigh"
  
  # 2/ (1) but rate K=20 resampling, rate L amp postfilter
  batch_process $fullfile "'amp_pf','phase_pf','K',20" "2_k20"

  # 3/ (1) with 2dB^2 random quant noise
  batch_process $fullfile "'amp_pf','phase_pf','K',20,'noise',2" "3_2dB"

  # 4/ (1) plus 1 stage VQ
  batch_process $fullfile "'amp_pf','phase_pf','K',20, \
  'vq1','../build_linux/train_k20_vq1.f32'" "4_vq1"
  
  # 5/ (1) plus 1 stage VQ with weighted search
  batch_process $fullfile "'amp_pf','phase_pf','K',20, \
  'vq1','../build_linux/train_k20_vq1.f32', \
  'weights', [ones(1,10) 0.5*ones(1,10)]" "5_vq1w"
   
  # 6/ (1) plus 2 stage VQ with weighted search
  batch_process $fullfile "'amp_pf','phase_pf','K',20, \
  'vq1','../build_linux/train_k20_vq1.f32', \
  'vq2','../build_linux/train_k20_vq2.f32', \
  'weights', [ones(1,10) 0.5*ones(1,10)]" "6_vq2w"
  
  # Codec 2 3200 & 700C controls
  c2enc 3200 $fullfile - | c2dec 3200 - - | sox -t .s16 -r 8000 -c 1 - ${out_dir}/${filename}_7_3200.wav
  c2enc 700C $fullfile - | c2dec 700C - - | sox -t .s16 -r 8000 -c 1 - ${out_dir}/${filename}_8_700C.wav
}

# Listen to effect of various eq algorithms.  Goal is to reduce dynamic range of
# data VQ has to handle
function postfilter_eq_test() {
  fullfile=$1
  filename=$(basename -- "$fullfile")
  extension="${filename##*.}"
  filename="${filename%.*}"
  mkdir -p $out_dir

  c2sim $fullfile --hpf --modelout ${filename}_model.bin

  # Amps Nb filtered, phase0, amp and phase postfilters, rate K
  echo "ratek3_batch; ratek3_batch_tool(\"${filename}\", \
       'A_out',\"${filename}_a.f32\",'H_out',\"${filename}_h.f32\", \
       'amp_pf','phase_pf','rateK'); quit;" \
  | octave -p ${CODEC2_PATH}/octave -qf
  c2sim $fullfile --hpf --phase0 --postfilter --amread ${filename}_a.f32 --hmread ${filename}_h.f32 -o - | \
      sox -t .s16 -r 8000 -c 1 - ${out_dir}/${filename}_1_Nb_p0_pf_rateK.wav

  # Compress dynamic range in post filter version 1
  echo "ratek3_batch; ratek3_batch_tool(\"${filename}\", \
       'A_out',\"${filename}_a.f32\",'H_out',\"${filename}_h.f32\", \
       'amp_pf','phase_pf','rateK', 'eq',1); quit;" \
  | octave -p ${CODEC2_PATH}/octave -qf
  c2sim $fullfile --hpf --phase0 --postfilter --amread ${filename}_a.f32 --hmread ${filename}_h.f32 -o - | \
      sox -t .s16 -r 8000 -c 1 - ${out_dir}/${filename}_2_eq1.wav

  # Compress dynamic range in post filter version 2
  echo "ratek3_batch; ratek3_batch_tool(\"${filename}\", \
       'A_out',\"${filename}_a.f32\",'H_out',\"${filename}_h.f32\", \
       'amp_pf','phase_pf','rateK', 'eq',2); quit;" \
  | octave -p ${CODEC2_PATH}/octave -qf
  c2sim $fullfile --hpf --phase0 --postfilter --amread ${filename}_a.f32 --hmread ${filename}_h.f32 -o - | \
      sox -t .s16 -r 8000 -c 1 - ${out_dir}/${filename}_3_eq2.wav
      
  # Compress dynamic range of B directly
  echo "ratek3_batch; ratek3_batch_tool(\"${filename}\", \
       'A_out',\"${filename}_a.f32\",'H_out',\"${filename}_h.f32\", \
       'amp_pf','phase_pf','rateK', 'DR', 40); quit;" \
   | octave -p ${CODEC2_PATH}/octave -qf
 c2sim $fullfile --hpf --phase0 --postfilter --amread ${filename}_a.f32 --hmread ${filename}_h.f32 -o - | \
      sox -t .s16 -r 8000 -c 1 - ${out_dir}/${filename}_4_dr.wav
       
  # Codec 2 3200 control
  c2enc 3200 $fullfile - | c2dec 3200 - - | sox -t .s16 -r 8000 -c 1 - ${out_dir}/${filename}_5_3200.wav
}

# Test postfilter/Am resampling with rate Lhigh and rate K processing, using ratek3_batch processing tool
# usage:
#   cd ~/codec2/build_linux
#   ../script/ratek_resampler.sh postfilter_rate_test

function postfilter_rate_test() {
  fullfile=$1
  filename=$(basename -- "$fullfile")
  extension="${filename##*.}"
  filename="${filename%.*}"
  mkdir -p $out_dir

  # orig amp and phase
  c2sim $fullfile --hpf -o - | sox -t .s16 -r 8000 -c 1 - ${out_dir}/${filename}_1_out.wav
  
  # orig amp and phase0.  Note uses c2sim internal Am->Hm, rather than our Octave version, bypassing Nb filtering
  c2sim $fullfile --hpf --phase0 --postfilter --modelout ${filename}_model.bin -o - | sox -t .s16 -r 8000 -c 1 - ${out_dir}/${filename}_2_p0.wav

  # amps Nb filtered, original phase
  echo "ratek3_batch; ratek3_batch_tool(\"${filename}\",'A_out',\"${filename}_a.f32\"); quit;" \
  | octave -p ${CODEC2_PATH}/octave -qf
  c2sim $fullfile --hpf --amread ${filename}_a.f32 -o - | \
      sox -t .s16 -r 8000 -c 1 - ${out_dir}/${filename}_3_Nb.wav

  # amps Nb filtered, phase0
  echo "ratek3_batch; ratek3_batch_tool(\"${filename}\",'A_out',\"${filename}_a.f32\",'H_out',\"${filename}_h.f32\"); quit;" \
  | octave -p ${CODEC2_PATH}/octave -qf
  c2sim $fullfile --hpf --phase0 --postfilter --amread ${filename}_a.f32 --hmread ${filename}_h.f32 -o - | \
      sox -t .s16 -r 8000 -c 1 - ${out_dir}/${filename}_4_Nb_p0.wav

  # Amps Nb filtered, phase0, amp and phase postfilters, rate Lhigh
  echo "ratek3_batch; ratek3_batch_tool(\"${filename}\",'A_out',\"${filename}_a.f32\",'H_out',\"${filename}_h.f32\",'amp_pf','phase_pf'); quit;" \
  | octave -p ${CODEC2_PATH}/octave -qf
  c2sim $fullfile --hpf --phase0 --postfilter --amread ${filename}_a.f32 --hmread ${filename}_h.f32 -o - | \
      sox -t .s16 -r 8000 -c 1 - ${out_dir}/${filename}_5_Nb_p0_pf_rateLhigh.wav

  # Amps Nb filtered, phase0, amp and phase postfilters, rate K
  echo "ratek3_batch; ratek3_batch_tool(\"${filename}\",'A_out',\"${filename}_a.f32\",'H_out',\"${filename}_h.f32\",'amp_pf','phase_pf','rateK'); quit;" \
  | octave -p ${CODEC2_PATH}/octave -qf
  c2sim $fullfile --hpf --phase0 --postfilter --amread ${filename}_a.f32 --hmread ${filename}_h.f32 -o - | \
      sox -t .s16 -r 8000 -c 1 - ${out_dir}/${filename}_6_Nb_p0_pf_rateK.wav

  # Codec 2 3200 control
  c2enc 3200 $fullfile - | c2dec 3200 - - | sox -t .s16 -r 8000 -c 1 - ${out_dir}/${filename}_7_3200.wav
}

# Process sample with various postfilter methods
# usage:
#   cd ~/codec2/build_linux
#   ../script/ratek_resampler.sh
function postfilter_test() {
  fullfile=$1
  filename=$(basename -- "$fullfile")
  extension="${filename##*.}"
  filename="${filename%.*}"
  mkdir -p $out_dir

  c2sim $fullfile --hpf -o - | sox -t .s16 -r 8000 -c 1 - ${out_dir}/${filename}_1_out.wav
  # TODO: uses c2sim internal Am->Hm, rather than our Octave version bypassing filtering
  c2sim $fullfile --hpf --phase0 --postfilter --dump $filename -o - | sox -t .s16 -r 8000 -c 1 - ${out_dir}/${filename}_2_p0.wav

  echo "ratek2_batch; ratek2_model_postfilter(\"${filename}\",\"${filename}_am.f32\"); quit;" \
  | octave -p ${CODEC2_PATH}/octave -qf
  c2sim $fullfile --hpf --amread ${filename}_am.f32 -o - | \
      sox -t .s16 -r 8000 -c 1 - ${out_dir}/${filename}_3_ratek.wav

  echo "ratek2_batch; ratek2_model_postfilter(\"${filename}\",\"${filename}_am.f32\",\"${filename}_hm.f32\"); quit;" \
  | octave -p ${CODEC2_PATH}/octave -qf
  c2sim $fullfile --hpf --phase0 --postfilter --amread ${filename}_am.f32 --hmread ${filename}_hm.f32 -o - | \
      sox -t .s16 -r 8000 -c 1 - ${out_dir}/${filename}_4_ratek_p0.wav

  echo "ratek2_batch; ratek2_model_postfilter(\"${filename}\",\"${filename}_am.f32\",\"${filename}_hm.f32\",1,0); quit;" \
  | octave -p ${CODEC2_PATH}/octave -qf
  c2sim $fullfile --hpf --phase0 --postfilter --amread ${filename}_am.f32 --hmread ${filename}_hm.f32 -o - | \
      sox -t .s16 -r 8000 -c 1 - ${out_dir}/${filename}_5_ratek_pf_p0.wav

  echo "ratek2_batch; ratek2_model_postfilter(\"${filename}\",\"${filename}_am.f32\",\"${filename}_hm.f32\",0,1); quit;" \
  | octave -p ${CODEC2_PATH}/octave -qf
  c2sim $fullfile --hpf --phase0 --postfilter --amread ${filename}_am.f32 --hmread ${filename}_hm.f32 -o - | \
      sox -t .s16 -r 8000 -c 1 - ${out_dir}/${filename}_6_ratek_p0_pf.wav

  echo "ratek2_batch; ratek2_model_postfilter(\"${filename}\",\"${filename}_am.f32\",\"${filename}_hm.f32\",1,1); quit;" \
  | octave -p ${CODEC2_PATH}/octave -qf
  c2sim $fullfile --hpf --phase0 --postfilter --amread ${filename}_am.f32 --hmread ${filename}_hm.f32 -o - | \
      sox -t .s16 -r 8000 -c 1 - ${out_dir}/${filename}_7_ratek_pf_p0_pf.wav

  echo "ratek2_batch; ratek2_model_postfilter(\"${filename}\",\"${filename}_am.f32\",\"${filename}_hm.f32\",0,0,1); quit;" \
  | octave -p ${CODEC2_PATH}/octave -qf
  c2sim $fullfile --hpf --phase0 --postfilter --amread ${filename}_am.f32 --hmread ${filename}_hm.f32 -o - | \
      sox -t .s16 -r 8000 -c 1 - ${out_dir}/${filename}_8_ratek_p0_smear.wav

  c2enc 3200 $fullfile - | c2dec 3200 - - | sox -t .s16 -r 8000 -c 1 - ${out_dir}/${filename}_9_3200.wav
}


function vq_test_subset() {
  fullfile=$1
  filename=$(basename -- "$fullfile")
  extension="${filename##*.}"
  filename="${filename%.*}"
  mkdir -p $out_dir
  vq1="train_b_vq1.f32"
  vq2="train_b_vq2.f32"
  vq1eq="train_eq_b_vq1.f32"
  vq2eq="train_eq_b_vq2.f32"

  c2sim $fullfile --hpf --modelout ${filename}_model.bin
  
  # Amps Nb filtered, phase0, amp and phase postfilters, rate K
  echo "ratek3_batch; ratek3_batch_tool(\"${filename}\",'A_out',\"${filename}_a.f32\",'H_out',\"${filename}_h.f32\",'amp_pf','phase_pf','rateK'); quit;" \
  | octave -p ${CODEC2_PATH}/octave -qf
  c2sim $fullfile --hpf --phase0 --postfilter --amread ${filename}_a.f32 --hmread ${filename}_h.f32 -o - | \
      sox -t .s16 -r 8000 -c 1 - ${out_dir}/${filename}_1_Nb_p0_pf_rateK.wav

  # As above plus subset
  echo "ratek3_batch; ratek3_batch_tool(\"${filename}\",\
       'A_out',\"${filename}_a.f32\",'H_out',\"${filename}_h.f32\",'amp_pf','phase_pf','rateK',
       'subset'); quit;" \
  | octave -p ${CODEC2_PATH}/octave -qf
  c2sim $fullfile --hpf --phase0 --postfilter --amread ${filename}_a.f32 --hmread ${filename}_h.f32 -o - | \
      sox -t .s16 -r 8000 -c 1 - ${out_dir}/${filename}_2_subset.wav

  # As above plus stage1 VQ
  echo "ratek3_batch; ratek3_batch_tool(\"${filename}\", \
        'A_out',\"${filename}_a.f32\",'H_out',\"${filename}_h.f32\",'amp_pf','phase_pf','rateK', \
        'vq1', \"${vq1}\",'subset'); quit;" \
        | octave -p ${CODEC2_PATH}/octave -qf
  c2sim $fullfile --hpf --phase0 --postfilter --amread ${filename}_a.f32 --hmread ${filename}_h.f32 -o - | \
      sox -t .s16 -r 8000 -c 1 - ${out_dir}/${filename}_3_vq1.wav

  # As above plus stage2 VQ
  echo "ratek3_batch; ratek3_batch_tool(\"${filename}\", \
        'A_out',\"${dr}${filename}_a.f32\",'H_out',\"${filename}_h.f32\",'amp_pf','phase_pf','rateK', \
        'vq1', \"${vq1}\", 'vq2', \"${vq2}\",'subset'); quit;" \
        | octave -p ${CODEC2_PATH}/octave -qf
  c2sim $fullfile --hpf --phase0 --postfilter --amread ${filename}_a.f32 --hmread ${filename}_h.f32 -o - | \
      sox -t .s16 -r 8000 -c 1 - ${out_dir}/${filename}_4_vq2.wav
      
  # Amps Nb filtered, phase0, amp and phase postfilters, rate K, EQ1
  echo "ratek3_batch; ratek3_batch_tool(\"${filename}\",'A_out',\"${filename}_a.f32\",'H_out',\"${filename}_h.f32\", \
       'amp_pf','phase_pf','rateK','eq',1); quit;" \
  | octave -p ${CODEC2_PATH}/octave -qf
  c2sim $fullfile --hpf --phase0 --postfilter --amread ${filename}_a.f32 --hmread ${filename}_h.f32 -o - | \
      sox -t .s16 -r 8000 -c 1 - ${out_dir}/${filename}_5_eq1.wav

  # As above plus subset
  echo "ratek3_batch; ratek3_batch_tool(\"${filename}\",'A_out',\"${filename}_a.f32\",'H_out',\"${filename}_h.f32\", \
       'amp_pf','phase_pf','rateK','eq',1,'subset'); quit;" \
  | octave -p ${CODEC2_PATH}/octave -qf
  c2sim $fullfile --hpf --phase0 --postfilter --amread ${filename}_a.f32 --hmread ${filename}_h.f32 -o - | \
      sox -t .s16 -r 8000 -c 1 - ${out_dir}/${filename}_6_eq1_subset.wav
  
  # As above plus stage1 VQ
  echo "ratek3_batch; ratek3_batch_tool(\"${filename}\", \
        'A_out',\"${filename}_a.f32\",'H_out',\"${filename}_h.f32\",'amp_pf','phase_pf','rateK', \
        'vq1', \"${vq1eq}\",'subset','eq',1); quit;" \
        | octave -p ${CODEC2_PATH}/octave -qf
  c2sim $fullfile --hpf --phase0 --postfilter --amread ${filename}_a.f32 --hmread ${filename}_h.f32 -o - | \
      sox -t .s16 -r 8000 -c 1 - ${out_dir}/${filename}_7_eq1_vq1.wav

  # As above plus stage2 VQ
  echo "ratek3_batch; ratek3_batch_tool(\"${filename}\", \
        'A_out',\"${filename}_a.f32\",'H_out',\"${filename}_h.f32\",'amp_pf','phase_pf','rateK', \
        'vq1', \"${vq1eq}\", 'vq2', \"${vq2eq}\",'subset','eq',1); quit;" \
        | octave -p ${CODEC2_PATH}/octave -qf
  c2sim $fullfile --hpf --phase0 --postfilter --amread ${filename}_a.f32 --hmread ${filename}_h.f32 -o - | \
      sox -t .s16 -r 8000 -c 1 - ${out_dir}/${filename}_8_eq1_vq2.wav
   
  # Codec 2 3200 & 700C controls
  c2enc 3200 $fullfile - | c2dec 3200 - - | sox -t .s16 -r 8000 -c 1 - ${out_dir}/${filename}_9_3200.wav
  c2enc 700C $fullfile - | c2dec 700C - - | sox -t .s16 -r 8000 -c 1 - ${out_dir}/${filename}_10_700C.wav
}

# generate amp postfiltered rate K training material from source speech file 
function gen_train() {
  fullfile=$1
  filename=$(basename -- "$fullfile")
  extension="${filename##*.}"
  filename="${filename%.*}"
  
  filename_b=${filename}_b.f32
  if [ $# -eq 2 ]; then
    filename_b=$2
  fi

  c2sim $fullfile --hpf --modelout ${filename}_model.bin
  echo "ratek3_batch; ratek3_batch_tool(\"${filename}\",'B_out',\"${filename_b}\", \
        ${options} \
        'K',${K}); quit;" \
  | octave -p ${CODEC2_PATH}/octave -qf
}

function gen_train_comp() {
  fullfile=$1
  filename=$(basename -- "$fullfile")
  extension="${filename##*.}"
  filename="${filename%.*}"
  
  filename_b=${filename}_b.f32
  if [ $# -eq 2 ]; then
    filename_b=$2
  fi

  c2sim $fullfile --hpf --prede --comp_gain 6 --comp 10000 --comp_dr 30 --modelout ${filename}_model.bin
  echo "ratek3_batch; ratek3_batch_tool(\"${filename}\",'B_out',\"${filename_b}\", \
        ${options} \
        'K',${K}); quit;" \
  | octave -p ${CODEC2_PATH}/octave -qf
}

function train_kmeans() {
  fullfile=$1
  filename=$(basename -- "$fullfile")
  extension="${filename##*.}"
  filename="${filename%.*}"
  res1=$(mktemp)
  
  # remove mean, train 2 stages - kmeans
  extract -t $K -s $Kst -e $Ken --lower 10 $removemean --writeall $fullfile ${filename}_nomean.f32
  vqtrain ${filename}_nomean.f32 $K $M --st $Kst --en $Ken -s 1e-3 ${filename}_vq1.f32 -r ${res1} --used ${filename}_used1.txt > kmeans_res1.txt
  vqtrain ${res1} $K $M --st $Kst --en $Ken  -s 1e-3 ${filename}_vq2.f32 -r res2.f32 --used ${filename}_used2.txt > kmeans_res2.txt
#  cat ${filename}_nomean.f32 | vq_mbest --mbest 5 -k $K --st $Kst --en $Ken  -q ${filename}_vq1.f32,${filename}_vq2.f32 >> /dev/null
}

# Generate a VQ by simply sampling ttrain database, idea is it avoids averaging
function train_no() {
  fullfile=$1
  filename=$(basename -- "$fullfile")
  extension="${filename##*.}"
  filename="${filename%.*}"
  res1=$(mktemp)
  
  filename_out=${filename}_no
  if [ $# -eq 2 ]; then
    filename_out=$2
  fi

  # remove mean, train 1 stages 
  extract -t $K -s $Kst -e $Ken --lower 10 --removemean --writeall $fullfile ${filename}_nomean.f32
  vqtrain ${filename}_nomean.f32 $K $M --st $Kst --en $Ken -s 1e-3 --notrain ${filename_out}_vq1.f32 -r ${res1} > ${filename_out}_res1.txt
  vqtrain ${res1} $K $M --st $Kst --en $Ken  -s 1e-3 --notrain ${filename_out}_vq2.f32 > ${filename_out}_res2.txt
  cat ${filename}_nomean.f32 | vq_mbest --mbest 5 -k $K --st $Kst --en $Ken  -q ${filename_out}_vq1.f32,${filename_out}_vq2.f32 >> /dev/null
}

function train_test() {
  fullfile=$1
  filename=$(basename -- "$fullfile")
  extension="${filename##*.}"
  filename="${filename%.*}"
  res1=$(mktemp)
  
  cat ${filename}_nomean.f32 | vq_mbest --mbest 5 -k $K --st $Kst --en $Ken -q $2 >> /dev/null        
}
        
function log2 {
    local x=0
    for (( y=$1-1 ; $y > 0; y >>= 1 )) ; do
        let x=$x+1
    done
    echo $x
}

function train_lbg() {
  fullfile=$1
  filename=$(basename -- "$fullfile")
  extension="${filename##*.}"
  filename="${filename%.*}"

  filename_out=${filename}_lbg
  if [ $# -eq 2 ]; then
    filename_out=$2
  fi
  
  # remove mean, extract columns from training data
  extract -t $K -s $Kst -e $Ken --lower $lower $removemean $meanl2 \
  --dynamicrange $dr $extract_options $drlate --writeall $fullfile ${filename_out}_nomean.f32

  # train 2 stages - LBG
  vqtrain ${filename_out}_nomean.f32 $K $M --st $Kst --en $Ken -s 1e-3 ${filename_out}_vq1.f32 -r res1.f32 --split > ${filename_out}_res1.txt
  if [ "$stage2" == "yes" ]; then
    vqtrain res1.f32 $K $M --st $Kst --en $Ken -s 1e-3 ${filename_out}_vq2.f32 -r res2.f32 --split > ${filename_out}_res2.txt
  fi
  if [ "$stage3" == "yes" ]; then
    vqtrain res2.f32 $K $M --st $Kst --en $Ken -s 1e-3 ${filename_out}_vq3.f32 --split > ${filename_out}_res3.txt
  fi
      
  # optionally compare stage3 search with mbest
  if [ "$mbest" == "yes" ]; then
    tmp=$(mktemp)
    results=${filename_out}_mbest3.txt
    rm ${results}
    log2M=$(log2 $M)
    for alog2M in $(seq 1 $log2M)
    do
      aM=$(( 2 ** $alog2M ))
      vqtrain res2.f32 $K $aM --st $Kst --en $Ken -s 1e-3 ${filename_out}_vq3.f32 --split > /dev/null
      cat ${filename_out}_nomean.f32 | \
          vq_mbest --mbest 5 -k $K -q ${filename_out}_vq1.f32,${filename_out}_vq2.f32,${filename_out}_vq3.f32 2>${tmp} >> /dev/null
      echo -n "$aM " >> ${results}
      cat ${tmp} | grep var | cut -d' ' -f 2 >> ${results}
    done
  fi

}

# Split VQ across freq |Kst..Ksp| and |Ksp+1 .. Ken|
function train_lbg_split() {
  fullfile=$1
  filename=$(basename -- "$fullfile")
  extension="${filename##*.}"
  filename="${filename%.*}"
  
  filename_out=${filename}_lbg
  if [ $# -eq 2 ]; then
    filename_out=$2
  fi
  
  # extract columns, remove mean and low energy vectors
  extract -t $K -s $Kst -e $Ken --lower 10 --removemean $fullfile ${filename_out}_nomean.f32
 
  # train 2 x 1 stage split VQ
  vqtrain ${filename_out}_nomean.f32 $K $M --st $Kst --en $Ksp -s 1e-3 ${filename_out}_vq1.f32 --split > ${filename_out}_res1.txt
  vqtrain ${filename_out}_nomean.f32 $K $M --st $(( $Ksp+1 )) --en $Ken -s 1e-3 ${filename_out}_vq2.f32 --split > ${filename_out}_res2.txt
}

# Split VQ across freq |Kst..Ksp| and |Ksp+1 .. Ken| and time (2 x 20ms)
function train_lbg_split_time() {
  fullfile=$1
  filename=$(basename -- "$fullfile")
  extension="${filename##*.}"
  filename="${filename%.*}"
  tmp=$(mktemp)
  tmp1=$(mktemp)
  tmp2=$(mktemp)
 
  filename_out=${filename}_lbg
  if [ $# -eq 2 ]; then
    filename_out=$2
  fi
  
  # extract all columns used, at 20ms time steps
  extract -t $K -s $Kst -e $Ken --timestep 2 $fullfile $tmp
  # Reshaping to vectors 2K long that include freq and 40ms of time, remove low
  # energy vectors across entire 40ms block
  extract -t $(( 2*K )) --lower 10 $tmp $tmp1
  # Reshaping back to vectors K long, extract mean across freq, implying
  # two means/40ms block that need quantising as side information
  extract -t $K --removemean $tmp1 $tmp2
  # Now extract two splits
  extract -t $K -s $Kst -e $Ksp $tmp2 ${filename_out}_nomean1.f32
  extract -t $K -s $(( $Ksp+1 )) -e $Ken $tmp2 ${filename_out}_nomean2.f32
  
  # train 2 x 1 stage split VQ, over 40ms time blocks
  K1=$(( $Ksp-$Kst+1 ))
  K2=$(( $Ken-$Ksp ))
  vqtrain ${filename_out}_nomean1.f32 $(( 2*K1 )) $M -s 1e-3 ${filename_out}_vq1.f32 --split > ${filename_out}_res1.txt
  vqtrain ${filename_out}_nomean2.f32 $(( 2*K2 )) $M -s 1e-3 ${filename_out}_vq2.f32 --split > ${filename_out}_res2.txt
 }

# predictive, 20ms steps
function train_lbg_pred() {
  fullfile=$1
  filename=$(basename -- "$fullfile")
  filename="${filename%.*}"
  tmp=$(mktemp)
  tmp1=$(mktemp)
  tmp2=$(mktemp)
 
  filename_out=${filename}_lbg
  if [ $# -eq 2 ]; then
    filename_out=$2
  fi
  
  # extract all columns used, at 20ms time steps
  extract -t $K -s $Kst -e $Ken --timestep 2 $fullfile $tmp
  # Reshaping to vectors 2K long that are 40ms of time, remove low
  # energy vectors across entire 40ms block
  extract -t $(( 2*K )) --lower 10 $tmp $tmp1
  # Reshaping back to vectors K long, and do prediction on 20ms steps
  extract -t $K -p 0.8 $tmp1 ${filename_out}_pred.f32
  
  # train 2 stage VQ, over 20ms time blocks
  vqtrain ${filename_out}_pred.f32 $K $M -s 1e-3 ${filename_out}_vq1.f32 -r res1.f32 --split > ${filename_out}_res1.txt
  vqtrain res1.f32 $K $M -s 1e-3 ${filename_out}_vq2.f32 --split > ${filename_out}_res2.txt
 }

# Try training with two different Nb
function train_Nb() {
  fullfile1=$1
  filename1=$(basename -- "$fullfile1")
  extension1="${filename1##*.}"
  filename1="${filename1%.*}"

  fullfile2=$2
  filename2=$(basename -- "$fullfile2")
  extension2="${filename2##*.}"
  filename2="${filename2%.*}"

  Nb1=$3
  Nb2=$4

  # remove mean, train 2 stages - LBG
  extract -t $K -s $Kst -e $Ken --removemean --writeall $fullfile1 ${filename1}_nomean.f32
  vqtrain ${filename1}_nomean.f32 $K $M  --st $Kst --en $Ken -s 1e-3 vq_stage1.f32 -r res1.f32 --split > lbg_res1.txt
  vqtrain res1.f32 $K $M  --st $Kst --en $Ken  -s 1e-3 vq_stage2.f32 -r res2.f32 --split > lbg_res2.txt

  extract -t $K -s $Kst -e $Ken --removemean --writeall $fullfile2 ${filename2}_nomean.f32
  vqtrain ${filename2}_nomean.f32 $K $M  --st $Kst --en $Ken -s 1e-3 vq_stage1.f32 -r res1.f32 --split > lbg_res3.txt
  vqtrain res1.f32 $K $M  --st $Kst --en $Ken  -s 1e-3 vq_stage2.f32 -r res2.f32 --split > lbg_res4.txt

  echo "lbg1=load('lbg_res1.txt'); lbg2=load('lbg_res2.txt'); \
        lbg3=load('lbg_res3.txt'); lbg4=load('lbg_res4.txt'); \
        hold on; \
        plot(log2(lbg1(:,1)),lbg1(:,2),'+-'); plot(log2(lbg2(:,1)),lbg2(:,2),'+-'); \
        plot(log2(lbg3(:,1)),lbg3(:,2),'o-'); plot(log2(lbg4(:,1)),lbg4(:,2),'o-'); \
        hold off; \
        leg = {'Nb=${Nb1} stage1','Nb=${Nb1} stage2','Nb=${Nb2} stage1','Nb=${Nb2} stage2'}; \
        h = legend(leg); legend('boxoff'); \
        set(gca, 'FontSize', 16); set (h, 'fontsize', 16);
        xlabel('Bits'); ylabel('Eq dB^2'); grid; \
        print(\"${Nb1}_${Nb2}_vq.png\",'-dpng','-S500,500'); \
        quit" | octave  -qf
}

if [ $# -gt 0 ]; then
  case $1 in
    postfilter_rate_test)
        postfilter_rate_test ../raw/big_dog.raw
        postfilter_rate_test ../raw/hts1a.raw
        postfilter_rate_test ../raw/hts2a.raw
        postfilter_rate_test ../raw/two_lines.raw
        postfilter_rate_test ../raw/cq_ref.raw
        postfilter_rate_test ../raw/kristoff.raw
        postfilter_rate_test ../raw/mmt1.raw      
        ;;
    gen_train)
        gen_train $2 $3
        ;;
    gen_train_comp)
        gen_train_comp $2 $3
        ;;
    train_kmeans)
        train_kmeans $2
        ;;
    train_test)
        train_test $2 $3
        ;;
    train_lbg)
        train_lbg $2 $3
        ;;
    train_lbg_split)
        train_lbg_split $2 $3
        ;;
    train_lbg_split_time)
        train_lbg_split_time $2 $3
        ;;
    train_lbg_pred)
        train_lbg_pred $2 $3
        ;;
    train_no)
        train_no $2 $3
        ;;
    vq_test)
        vq_test ../raw/big_dog.raw
        vq_test ../raw/two_lines.raw
        ;;
    vq_test_230202)
        vq_test_230202 ../raw/big_dog.raw
        vq_test_230202 ../raw/two_lines.raw
        ;;
    dr_vq_test_230204)
        dr_vq_test_230204 ../raw/big_dog.raw
        dr_vq_test_230204 ../raw/two_lines.raw
        ;;
    vq_test_230204)
        #rm -f ${out_dir}/zlog.txt
        #vq_test_230204 ../raw/forig.raw
        #vq_test_230204 ../raw/big_dog.raw
        #vq_test_230204 ../raw/two_lines.raw
        #vq_test_230204 ../raw/cq_ref.raw
        #vq_test_230204 ../raw/morig.raw
        #vq_test_230204 ../raw/hts2a.raw        
        #vq_test_230204  ../wav/vk5dgr_testing_8k.wav  
        vq_test_230204  ../raw/ship.raw
        vq_test_230204  ../raw/sickness.raw
        vq_test_230204  ../raw/pickle.raw
        vq_test_230204  ../raw/tea.raw
        ;;
        
    vq_test_230217)
        #rm -f ${out_dir}/zlog.txt
        vq_test_230217 ../raw/morig.raw
        #vq_test_230217 ../raw/big_dog.raw
        vq_test_230217 ../raw/cq_ref.raw
        vq_test_230217 ../raw/two_lines.raw
        #vq_test_230217 ../raw/hts1a.raw        
        #vq_test_230217 ../raw/hts2a.raw        
        vq_test_230217 ../raw/kristoff.raw        
        vq_test_230217 ../raw/forig.raw     
        vq_test_230217 ../raw/ve9qrp_10s.raw     
        ;;

    vq_test_230226)
        vq_test_230326 ../raw/forig.raw     
        vq_test_230226 ../raw/big_dog.raw
        vq_test_230226 ../raw/two_lines.raw
        vq_test_230226 ../raw/pickle.raw
        vq_test_230226 ../raw/tea.raw
        ;;
 
    comp_test_230323)
        comp_test_230323 ../raw/big_dog.raw
        comp_test_230323 ../raw/forig.raw     
        comp_test_230323 ../raw/two_lines.raw
        comp_test_230323 ../raw/pickle.raw
        comp_test_230323 ../raw/tea.raw
        comp_test_230323 ../raw/kristoff.raw        
        comp_test_230323 ../raw/ve9qrp_10s.raw     
        comp_test_230323 ../raw/mmt1.raw     
        ;;
 
    vq_test_230331)
        vq_test_230331 ../raw/big_dog.raw
        #vq_test_230331 ../raw/forig.raw     
        #vq_test_230331 ../raw/two_lines.raw
        #vq_test_230331 ../raw/pickle.raw
        #vq_test_230331 ../raw/tea.raw
        #vq_test_230331  ../raw/ship.raw
        #vq_test_230331  ../raw/sickness.raw
        #vq_test_230331 ../raw/kristoff.raw        
        #vq_test_230331 ../raw/ve9qrp_10s.raw     
        #vq_test_230331 ../raw/mmt1.raw     
        ;;
        
    vq_test_subset)
        vq_test_subset ../raw/big_dog.raw
        vq_test_subset ../raw/two_lines.raw
        ;;
    postfilter_eq_test)
        postfilter_eq_test ../raw/big_dog.raw
        postfilter_eq_test ../raw/two_lines.raw
        ;;
  esac
else
  echo "usage:
  echo "  ./ratek_resampler.sh command [options ...]""
fi

