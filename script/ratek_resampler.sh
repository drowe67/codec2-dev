#!/bin/bash -x
# ratek_resampler.sh
# David Rowe Sep 2022
#
# Support for rate K resampler experiments see doc/ratek_resampler

CODEC2_PATH=$HOME/codec2
PATH=$PATH:$CODEC2_PATH/build_linux/src:$CODEC2_PATH/build_linux/misc
K="${K:-30}"
M="${M:-4096}"
Kst="${Kst:-0}"
Ken="${Ken:-29}"
out_dir="${out_dir:-ratek_out}"
Nb=20

# Listen to effect of not restoring slope, which has potential to lower dynamic range
function postfilter_slope_test() {
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
      sox -t .s16 -r 8000 -c 1 - ${out_dir}/${filename}_1_slope.wav

  # As above but don't restore slope after amplitude post filtering
  echo "ratek3_batch; ratek3_batch_tool(\"${filename}\", \
       'A_out',\"${filename}_a.f32\",'H_out',\"${filename}_h.f32\", \
       'amp_pf','phase_pf','rateK', 'eq'); quit;" \
  | octave -p ${CODEC2_PATH}/octave -qf
  c2sim $fullfile --hpf --phase0 --postfilter --amread ${filename}_a.f32 --hmread ${filename}_h.f32 -o - | \
      sox -t .s16 -r 8000 -c 1 - ${out_dir}/${filename}_2_eq.wav

  # Codec 2 3200 control
  c2enc 3200 $fullfile - | c2dec 3200 - - | sox -t .s16 -r 8000 -c 1 - ${out_dir}/${filename}_3_3200.wav
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

# Process sample with various methods including 1 and 2 stage VQ
function vq_test() {
  fullfile=$1
  filename=$(basename -- "$fullfile")
  extension="${filename##*.}"
  filename="${filename%.*}"
  mkdir -p $out_dir

  c2sim $fullfile --hpf --modelout ${filename}_model.bin
  
  # Amps Nb filtered, phase0, amp and phase postfilters, rate K
  echo "ratek3_batch; ratek3_batch_tool(\"${filename}\",'A_out',\"${filename}_a.f32\",'H_out',\"${filename}_h.f32\",'amp_pf','phase_pf','rateK'); quit;" \
  | octave -p ${CODEC2_PATH}/octave -qf
  c2sim $fullfile --hpf --phase0 --postfilter --amread ${filename}_a.f32 --hmread ${filename}_h.f32 -o - | \
      sox -t .s16 -r 8000 -c 1 - ${out_dir}/${filename}_1_Nb_p0_pf_rateK.wav

  # As above plus stage1 VQ
  echo "ratek3_batch; ratek3_batch_tool(\"${filename}\", \
        'A_out',\"${filename}_a.f32\",'H_out',\"${filename}_h.f32\",'amp_pf','phase_pf','rateK', \
        'vq_stage1', 'vq_stage1.f32'); quit;" \
        | octave -p ${CODEC2_PATH}/octave -qf
  c2sim $fullfile --hpf --phase0 --postfilter --amread ${filename}_a.f32 --hmread ${filename}_h.f32 -o - | \
      sox -t .s16 -r 8000 -c 1 - ${out_dir}/${filename}_2_vq1.wav

  # As above plus stage2 VQ
  echo "ratek3_batch; ratek3_batch_tool(\"${filename}\", \
        'A_out',\"${filename}_a.f32\",'H_out',\"${filename}_h.f32\",'amp_pf','phase_pf','rateK', \
        'vq_stage1', 'vq_stage1.f32', 'vq_stage2', 'vq_stage2.f32'); quit;" \
        | octave -p ${CODEC2_PATH}/octave -qf
  c2sim $fullfile --hpf --phase0 --postfilter --amread ${filename}_a.f32 --hmread ${filename}_h.f32 -o - | \
      sox -t .s16 -r 8000 -c 1 - ${out_dir}/${filename}_3_vq2.wav

  # Codec 2 3200 & 700C controls
  c2enc 3200 $fullfile - | c2dec 3200 - - | sox -t .s16 -r 8000 -c 1 - ${out_dir}/${filename}_4_3200.wav
  c2enc 700C $fullfile - | c2dec 700C - - | sox -t .s16 -r 8000 -c 1 - ${out_dir}/${filename}_5_700C.wav
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
  echo "ratek3_batch; ratek3_batch_tool(\"${filename}\",'B_out',\"${filename_b}\",'amp_pf','K',${K}); quit;" \
  | octave -p ${CODEC2_PATH}/octave -qf
}

function gen_train_eq() {
  fullfile=$1
  filename=$(basename -- "$fullfile")
  extension="${filename##*.}"
  filename="${filename%.*}"
  
  c2sim $fullfile --hpf --dump $filename
  echo "ratek3_batch; ratek3_batch_tool(\"${filename}\",'B_out',\"${filename}_eq_b.f32\",'amp_pf','eq'); quit;" \
  | octave -p ${CODEC2_PATH}/octave -qf
}

function train_kmeans() {
  fullfile=$1
  filename=$(basename -- "$fullfile")
  extension="${filename##*.}"
  filename="${filename%.*}"
  res1=$(mktemp)
  
  # remove mean, train 2 stages - kmeans
  extract -t $K -s $Kst -e $Ken --lower 10 --removemean --writeall --dynamicrange 40 $fullfile ${filename}_nomean.f32
  vqtrain ${filename}_nomean.f32 $K $M --st $Kst --en $Ken -s 1e-3 ${filename}_vq1.f32 -r ${res1} --used ${filename}_used1.txt > kmeans_res1.txt
  vqtrain ${res1} $K $M --st $Kst --en $Ken  -s 1e-3 ${filename}_vq2.f32 -r res2.f32 --used ${filename}_used2.txt > kmeans_res2.txt
}

# comparing kmeans to lbg

function train_kmeans_lbg() {
  fullfile=$1
  filename=$(basename -- "$fullfile")
  extension="${filename##*.}"
  filename="${filename%.*}"

  # remove mean, train 2 stages - kmeans
  extract -t $K -s $Kst -e $Ken --removemean --writeall $fullfile ${filename}_nomean.f32
  vqtrain ${filename}_nomean.f32 $K $M  --st $Kst --en $Ken -s 1e-3 vq_stage1.f32 -r res1.f32 > kmeans_res1.txt
  vqtrain res1.f32 $K $M  --st $Kst --en $Ken  -s 1e-3 vq_stage2.f32 -r res2.f32 > kmeans_res2.txt

  # remove mean, train 2 stages - LBG
  extract -t $K -s $Kst -e $Ken --removemean --writeall $fullfile ${filename}_nomean.f32
  vqtrain ${filename}_nomean.f32 $K $M  --st $Kst --en $Ken -s 1e-3 vq_stage1.f32 -r res1.f32 --split > lbg_res1.txt
  vqtrain res1.f32 $K $M  --st $Kst --en $Ken  -s 1e-3 vq_stage2.f32 -r res2.f32 --split > lbg_res2.txt
  cat ${filename}_nomean.f32 | vq_mbest --mbest 5 -k $K -q vq_stage1.f32,vq_stage2.f32 >> /dev/null

  echo "kmeans1=load('kmeans_res1.txt'); kmeans2=load('kmeans_res2.txt'); \
        lbg1=load('lbg_res1.txt'); lbg2=load('lbg_res2.txt'); \
        hold on; \
        plot(log2(kmeans1(:,1)),kmeans1(:,2),'+-','markersize', 15); plot(log2(kmeans2(:,1)),kmeans2(:,2),'+-','markersize', 15); \
        plot(log2(lbg1(:,1)),lbg1(:,2),'+-'); plot(log2(lbg2(:,1)),lbg2(:,2),'+-'); \
        hold off; \
        leg = {'kmeans stage1','kmeans stage2','lbg stage1','lbg stage2'}; \
        h = legend(leg); legend('boxoff'); \
        set(gca, 'FontSize', 16); set (h, 'fontsize', 16);
        xlabel('Bits'); ylabel('Eq dB^2'); grid; \
        print(\"${filename}_vq.png\",'-dpng','-S500,500'); \
        quit" | octave  -qf
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
    gen_train_eq)
        gen_train_eq $2
        ;;
    vq_train)
        train_kmeans $2
        ;;
    vq_test)
        vq_test ../raw/big_dog.raw
        vq_test ../raw/two_lines.raw
        ;;
    postfilter_slope_test)
        postfilter_slope_test ../raw/big_dog.raw
        postfilter_slope_test ../raw/two_lines.raw
        ;;
  esac
else
  echo "usage:
  echo "  ./ratek_resampler.sh command [options ...]""
fi

# TODO: make these selectable via CLI
#postfilter_test ../raw/hts1a.raw
#postfilter_test ../raw/two_lines.raw
#postfilter_test ../raw/cq_ref.raw
#postfilter_test ../raw/kristoff.raw
#postfilter_test ../raw/mmt1.raw

#test $1

#train_kmeans $1

#../script/ratek_resampler.sh ../octave/train_120_Nb20_K30.f32
#train_kmeans_lbg $1

# ../script/ratek_resampler.sh ../octave/train_120_Nb20_K30.f32 ../octave/train_120_Nb100_K30.f32 20 100
#train_Nb $1 $2 $3 $4
