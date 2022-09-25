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

#../script/ratek_resampler.sh ../octave/train_120_Nb20_K30.f32
train_kmeans_lbg $1

# ../script/ratek_resampler.sh ../octave/train_120_Nb20_K30.f32 ../octave/train_120_Nb100_K30.f32 20 100 
#train_Nb $1 $2 $3 $4
