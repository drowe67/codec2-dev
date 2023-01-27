# ratek_resampler.mk
#
# Makefile to build up curves for various VQ schemes
#
# usage:
#   cd codec2/build/linux
#   make -f ../script/ratek_resampler.mk

TRAIN  ?= train_120
M ?= 512

HELL  := /bin/bash
CODEC2 := $(HOME)/codec2
TRAIN_FULL := ~/Downloads/$(TRAIN).spc

PLOT_DATA := $(TRAIN)_lbg_res1.txt $(TRAIN)_lbg_res2.txt $(TRAIN)_lbg_mbest2.txt \
             $(TRAIN)_sub_res1.txt $(TRAIN)_sub_res2.txt \
             $(TRAIN)_k20_res1.txt $(TRAIN)_k20_res2.txt \
             $(TRAIN)_splt_res1.txt $(TRAIN)_splt_res2.txt \
	     $(TRAIN)_stf_res1.txt $(TRAIN)_stf_res2.txt

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
             \"$(TRAIN)_stf_res2.txt\",'bk--+;k20 stf2;' \
             ); quit" | octave-cli -p $(CODEC2)/octave --no-init-file

# (1) no amp PF before VQ, include 2nd stage mbest
$(TRAIN)_lbg_res1.txt $(TRAIN)_lbg_res2.txt $(TRAIN)_lbg_mbest2.txt: $(TRAIN)_b.f32
	mbest="yes" M=$(M) ../script/ratek_resampler.sh train_lbg $(TRAIN)_b.f32 $(TRAIN)_lbg

# (2) as above, subset
$(TRAIN)_sub_res1.txt $(TRAIN)_sub_res2.txt: $(TRAIN)_b.f32
	Kst=2 Ken=24 M=$(M) ../script/ratek_resampler.sh train_lbg $(TRAIN)_b.f32 $(TRAIN)_sub

# as per (1), K=20
$(TRAIN)_k20_res1.txt $(TRAIN)_k20_res2.txt: $(TRAIN)_b20.f32
	K=20 Kst=0 Ken=19 M=$(M) ../script/ratek_resampler.sh train_lbg $(TRAIN)_b20.f32 $(TRAIN)_k20

# K=20 split, energy removed first
$(TRAIN)_splt_res1.txt $(TRAIN)_splt_res2.txt: $(TRAIN)_b20.f32
	K=20 Kst=0 Ksp=11 Ken=19  M=$(M) ../script/ratek_resampler.sh train_lbg_split $(TRAIN)_b20.f32 $(TRAIN)_splt

# K=20 split, time and freq, energy removed first, slight subset at HF end
$(TRAIN)_stf_res1.txt $(TRAIN)_stf_res2.txt: $(TRAIN)_b20.f32
	K=20 Kst=0 Ksp=11 Ken=19 M=$(M) ../script/ratek_resampler.sh train_lbg_split_time $(TRAIN)_b20.f32 $(TRAIN)_stf

$(TRAIN)_b.f32:
	../script/ratek_resampler.sh gen_train $(TRAIN_FULL)

$(TRAIN)_b20.f32:
	K=20 ../script/ratek_resampler.sh gen_train $(TRAIN_FULL) $(TRAIN)_b20.f32

clean:
	rm -f $(PLOT_DATA)
	