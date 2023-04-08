# ratek_resampler.mk
#
# Makefile to build up curves for various VQ schemes
#
# usage:
#
# 1. Use train_120.spc (runs fairly quickly, good for initial tests)
#
#   cd codec2/build/linux
#   make -f ../script/ratek_resampler.mk
#
# 2. Use train.spc 
#
#    TRAIN=train M=4096 make -f ../script/ratek_resampler.mk
#
# 3. Generate .tex/.eps file for use in report
#
#  epslatex_path=../doc/ratek_resampler/ TRAIN=train M=4096 make -f ../script/ratek_resampler.mk

TRAIN ?= train_120
M ?= 512

SHELL  := /bin/bash
CODEC2 := $(HOME)/codec2
TRAIN_FULL := ~/Downloads/$(TRAIN).spc

PLOT_DATA := $(TRAIN)_lbg_res1.txt $(TRAIN)_lbg_res2.txt $(TRAIN)_lbg_mbest2.txt \
             $(TRAIN)_sub_res1.txt $(TRAIN)_sub_res2.txt \
             $(TRAIN)_k20_res1.txt $(TRAIN)_k20_res2.txt \
             $(TRAIN)_splt_res1.txt $(TRAIN)_splt_res2.txt \
	     $(TRAIN)_stf_res1.txt $(TRAIN)_stf_res2.txt

PLOT_DATA1 := $(TRAIN)_k20_res1.txt $(TRAIN)_k20_res2.txt \
	      $(TRAIN)_three_res1.txt $(TRAIN)_three_res2.txt \
	      $(TRAIN)_three_res3.txt $(TRAIN)_three_mbest3.txt
	      
all: $(TRAIN)_ratek.png $(TRAIN)_ratek1.png

# always re-render PNGs
FORCE:

$(TRAIN)_ratek.png: $(PLOT_DATA) FORCE
	echo "ratek_resampler_plot(\"$(TRAIN)_ratek.png\", \
             \"$(TRAIN)_lbg_res1.txt\",'g-+;lbg1;', \
             'continue',\"$(TRAIN)_lbg_res2.txt\",'g-+;lbg2;', \
             \"$(TRAIN)_sub_res1.txt\",'b-+;sub1;', \
             'continue',\"$(TRAIN)_sub_res2.txt\",'b-+;sub2;', \
             \"$(TRAIN)_k20_res1.txt\",'r-*;k20 1;', \
             'continue',\"$(TRAIN)_k20_res2.txt\",'r-*;k20 2;', \
             \"$(TRAIN)_splt_res1.txt\",'c-+;k20 splt 1;', \
             \"$(TRAIN)_splt_res2.txt\",'c--+;k20 splt 2;', \
             \"$(TRAIN)_stf_res1.txt\",'bk-+;k20 stf1;', \
             \"$(TRAIN)_stf_res2.txt\",'bk--+;k20 stf2;' \
             ); quit" | octave-cli -p $(CODEC2)/octave --no-init-file

$(TRAIN)_ratek1.png: $(PLOT_DATA1) FORCE
	echo "ratek_resampler_plot(\"$(TRAIN)_ratek1.png\", \
	     \"$(TRAIN)_k20_res1.txt\",'r-x;k20 12b 1;', \
             \"$(TRAIN)_three_res1.txt\",'b-+;k20 9b 1;', \
             'continue',\"$(TRAIN)_three_res2.txt\",'b-+;k20 9b 2;', \
             'continue',\"$(TRAIN)_three_res3.txt\",'b-+;k20 9b 3;' \
             ); quit" | octave-cli -p $(CODEC2)/octave --no-init-file

# (1) K=30 two stage VQ
$(TRAIN)_lbg_res1.txt $(TRAIN)_lbg_res2.txt $(TRAIN)_lbg_mbest2.txt: $(TRAIN)_b.f32
	mbest="yes" M=$(M) ../script/ratek_resampler.sh train_lbg $(TRAIN)_b.f32 $(TRAIN)_lbg

# (2) as above, subset
$(TRAIN)_sub_res1.txt $(TRAIN)_sub_res2.txt: $(TRAIN)_b.f32
	Kst=2 Ken=24 M=$(M) ../script/ratek_resampler.sh train_lbg $(TRAIN)_b.f32 $(TRAIN)_sub

# (3) K=20 two stage VQ
$(TRAIN)_k20_res1.txt $(TRAIN)_k20_res2.txt: $(TRAIN)_b20.f32
	K=20 Kst=0 Ken=19 M=$(M) ../script/ratek_resampler.sh train_lbg $(TRAIN)_b20.f32 $(TRAIN)_k20

# (4) K=20 split, energy removed first
$(TRAIN)_splt_res1.txt $(TRAIN)_splt_res2.txt: $(TRAIN)_b20.f32
	K=20 Kst=0 Ksp=11 Ken=19  M=$(M) ../script/ratek_resampler.sh train_lbg_split $(TRAIN)_b20.f32 $(TRAIN)_splt

# (5) K=20 split, time and freq, energy removed first
$(TRAIN)_stf_res1.txt $(TRAIN)_stf_res2.txt: $(TRAIN)_b20.f32
	K=20 Kst=0 Ksp=11 Ken=19 M=$(M) ../script/ratek_resampler.sh train_lbg_split_time $(TRAIN)_b20.f32 $(TRAIN)_stf

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

$(TRAIN)_b.f32:
	../script/ratek_resampler.sh gen_train $(TRAIN_FULL)

$(TRAIN)_b20.f32:
	K=20 ../script/ratek_resampler.sh gen_train $(TRAIN_FULL) $(TRAIN)_b20.f32

clean:
	rm -f $(PLOT_DATA)
	