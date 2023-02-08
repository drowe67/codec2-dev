# ratek_resampler.mk
#
# Makefile to build up curves for various VQ schemes
#
# usage:
#   cd codec2/build/linux
#   make -f ../script/ratek_resampler.mk

TRAIN  ?= train_120
M ?= 512

SHELL  := /bin/bash
CODEC2 := $(HOME)/codec2
TRAIN_FULL := ~/Downloads/$(TRAIN).spc

PLOT_DATA := $(TRAIN)_lbg_res1.txt $(TRAIN)_lbg_res2.txt $(TRAIN)_lbg_mbest2.txt \
             $(TRAIN)_sub_res1.txt $(TRAIN)_sub_res2.txt \
             $(TRAIN)_k20_res1.txt $(TRAIN)_k20_res2.txt \
             $(TRAIN)_splt_res1.txt $(TRAIN)_splt_res2.txt \
	     $(TRAIN)_stf_res1.txt $(TRAIN)_stf_res2.txt \
	     $(TRAIN)_comp_res1.txt $(TRAIN)_comp_res2.txt

PLOT_DATA1 := $(TRAIN)_k20_res1.txt $(TRAIN)_k20_res2.txt \
	      $(TRAIN)_pre_res1.txt $(TRAIN)_pre_res2.txt \
	      $(TRAIN)_comp_res1.txt $(TRAIN)_comp_res2.txt \
	      $(TRAIN)_three_res1.txt $(TRAIN)_three_res2.txt \
	      $(TRAIN)_three_res3.txt $(TRAIN)_three_mbest3.txt \
	      $(TRAIN)_pred_res1.txt $(TRAIN)_pred_res2.txt \
	      $(TRAIN)_sub_res1.txt $(TRAIN)_sub_res2.txt \
	      $(TRAIN)_suba_res1.txt $(TRAIN)_suba_res2.txt \
	      $(TRAIN)_pred2_res1.txt $(TRAIN)_pred2_res2.txt
	      
all: $(TRAIN)_ratek.png $(TRAIN)_ratek1.png

$(TRAIN)_ratek.png: $(PLOT_DATA)
	echo "ratek_resampler_plot(\"$(TRAIN)_ratek.png\", \
             \"$(TRAIN)_lbg_res1.txt\",'g-+;lbg1;', \
             'continue',\"$(TRAIN)_lbg_res2.txt\",'g-+;lbg2;', \
             'continue',\"$(TRAIN)_lbg_mbest2.txt\",'gx;mbest2;', \
             \"$(TRAIN)_sub_res1.txt\",'b-+;sub1;', \
             'continue',\"$(TRAIN)_sub_res2.txt\",'b-+;sub2;', \
             \"$(TRAIN)_k20_res1.txt\",'r-*;k20 1;', \
             'continue',\"$(TRAIN)_k20_res2.txt\",'r-*;k20 2;', \
             \"$(TRAIN)_splt_res1.txt\",'c-+;k20 splt 1;', \
             \"$(TRAIN)_splt_res2.txt\",'c--+;k20 splt 2;', \
             \"$(TRAIN)_stf_res1.txt\",'bk-+;k20 stf1;', \
             \"$(TRAIN)_stf_res2.txt\",'bk--+;k20 stf2;', \
             \"$(TRAIN)_comp_res1.txt\",'c-o;k20 comp1;', \
             'continue', \"$(TRAIN)_comp_res2.txt\",'c-o;k20 comp2;' \
             ); quit" | octave-cli -p $(CODEC2)/octave --no-init-file

$(TRAIN)_ratek1.png: $(PLOT_DATA1)
	echo "ratek_resampler_plot(\"$(TRAIN)_ratek1.png\", \
	     \"$(TRAIN)_k20_res1.txt\",'r-x;k20 1;', \
	     'continue',\"$(TRAIN)_k20_res2.txt\",'r-x;k20 2;', \
             \"$(TRAIN)_sub_res1.txt\",'b-+;sub1 2..24;', \
             'continue',\"$(TRAIN)_sub_res2.txt\",'b-+;sub2 2..24;', \
             \"$(TRAIN)_suba_res1.txt\",'g-+;suba1 0..17;', \
             'continue',\"$(TRAIN)_suba_res2.txt\",'g-+;suba2 0..17;', \
             \"$(TRAIN)_pred2_res1.txt\",'c-+;pred2 1;', \
             'continue',\"$(TRAIN)_pred2_res2.txt\",'c-+;pred2 2;' \
             ); quit" | octave-cli -p $(CODEC2)/octave --no-init-file

# (1) no amp PF before VQ, include 2nd stage mbest
$(TRAIN)_lbg_res1.txt $(TRAIN)_lbg_res2.txt $(TRAIN)_lbg_mbest2.txt: $(TRAIN)_b.f32
	mbest="yes" M=$(M) ../script/ratek_resampler.sh train_lbg $(TRAIN)_b.f32 $(TRAIN)_lbg

# (2) as above, subset
$(TRAIN)_sub_res1.txt $(TRAIN)_sub_res2.txt: $(TRAIN)_b.f32
	Kst=2 Ken=24 M=$(M) ../script/ratek_resampler.sh train_lbg $(TRAIN)_b.f32 $(TRAIN)_sub

# (3) K=20
$(TRAIN)_k20_res1.txt $(TRAIN)_k20_res2.txt: $(TRAIN)_b20.f32
	K=20 Kst=0 Ken=19 M=$(M) ../script/ratek_resampler.sh train_lbg $(TRAIN)_b20.f32 $(TRAIN)_k20

# (4) K=20 split, energy removed first
$(TRAIN)_splt_res1.txt $(TRAIN)_splt_res2.txt: $(TRAIN)_b20.f32
	K=20 Kst=0 Ksp=11 Ken=19  M=$(M) ../script/ratek_resampler.sh train_lbg_split $(TRAIN)_b20.f32 $(TRAIN)_splt

# (5) K=20 split, time and freq, energy removed first, slight subset at HF end
$(TRAIN)_stf_res1.txt $(TRAIN)_stf_res2.txt: $(TRAIN)_b20.f32
	K=20 Kst=0 Ksp=11 Ken=19 M=$(M) ../script/ratek_resampler.sh train_lbg_split_time $(TRAIN)_b20.f32 $(TRAIN)_stf

# (6) K=20, pre-emphasis on {Am} before rate K, VQ
$(TRAIN)_pre_res1.txt $(TRAIN)_pre_res2.txt: $(TRAIN)_b20_pre.f32
	K=20 Kst=0 Ken=19 M=$(M) ../script/ratek_resampler.sh train_lbg $(TRAIN)_b20_pre.f32 $(TRAIN)_pre

# (7) As per (6), but limit dynamic range
$(TRAIN)_comp_res1.txt $(TRAIN)_comp_res2.txt: $(TRAIN)_b20_pre.f32
	K=20 Kst=0 Ken=19 M=$(M) dr=30 ../script/ratek_resampler.sh train_lbg $(TRAIN)_b20_pre.f32 $(TRAIN)_comp

# (8) As per (7), but limit dynamic range after mean removal
$(TRAIN)_compl_res1.txt $(TRAIN)_compl_res2.txt: $(TRAIN)_b20_pre.f32
	K=20 Kst=0 Ken=19 M=$(M) dr=30 drlate="--drlate" ../script/ratek_resampler.sh train_lbg $(TRAIN)_b20_pre.f32 $(TRAIN)_compl

# (9) no training, just random sampling
$(TRAIN)_no_res1.txt $(TRAIN)_no_res2.txt: $(TRAIN)_b20.f32
	K=20 Kst=0 Ken=19 M=$(M) ../script/ratek_resampler.sh train_no $(TRAIN)_b20.f32 $(TRAIN)_no

# (10) three stages 9 bit/stage
$(TRAIN)_three_res1.txt $(TRAIN)_three_res2.txt $(TRAIN)_three_res3.txt $(TRAIN)_three_mbest3.txt: $(TRAIN)_b20.f32
	K=20 Kst=0 Ken=19 M=512 stage3="yes" mbest="yes" ../script/ratek_resampler.sh train_lbg $(TRAIN)_b20.f32 $(TRAIN)_three

# (11) predictive 20ms
$(TRAIN)_pred_res1.txt $(TRAIN)_pred_res2.txt: $(TRAIN)_b20.f32
	K=20 Kst=0 Ken=19 M=$(M) removemean=" " lower=-100 extract_options="--p 0.9 -d 2" \
	../script/ratek_resampler.sh train_lbg_pred $(TRAIN)_b20.f32 $(TRAIN)_pred

# (12) Another attempt at subset
$(TRAIN)_suba_res1.txt $(TRAIN)_suba_res2.txt: $(TRAIN)_b20.f32
	K=20 Kst=0 Ken=17 M=$(M) ../script/ratek_resampler.sh train_lbg $(TRAIN)_b20.f32 $(TRAIN)_suba

# (13) Another attempt at subset
$(TRAIN)_pred2_res1.txt $(TRAIN)_pred2_res2.txt: $(TRAIN)_b20.f32
	K=20 Kst=0 Ken=19 M=$(M) removemean=" " extract_options="-d 4 --pred2" \
	../script/ratek_resampler.sh train_lbg $(TRAIN)_b20.f32 $(TRAIN)_pred2


$(TRAIN)_b.f32:
	../script/ratek_resampler.sh gen_train $(TRAIN_FULL)

$(TRAIN)_b20.f32:
	K=20 ../script/ratek_resampler.sh gen_train $(TRAIN_FULL) $(TRAIN)_b20.f32

$(TRAIN)_b20_pre.f32:
	K=20 options="'pre'," ../script/ratek_resampler.sh gen_train $(TRAIN_FULL) $(TRAIN)_b20_pre.f32

clean:
	rm -f $(PLOT_DATA)
	