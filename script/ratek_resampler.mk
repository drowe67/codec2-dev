# ratek_resampler.mk
# usage:
#   cd codec2/build/linux
#   make -f ../script/ratek_resampler.mk

SHELL  := /bin/bash
CODEC2 := $(HOME)/codec2
TRAIN_FULL := ~/Downloads/train_120.spc
TRAIN      := train_120

ratek.png: $(TRAIN)_b_lbg_res1.txt $(TRAIN)_b_lbg_res2.txt $(TRAIN)_b_lbg_mbest2.txt
	echo "ratek_resampler_plot(\"ratek.png\", \
             \"$(TRAIN)_b_lbg_res1.txt\",'g-+;stage11;', \
             'continue',\"$(TRAIN)_b_lbg_res2.txt\",'r-+;stage2;', \
             'continue',\"$(TRAIN)_b_lbg_mbest2.txt\",'r-x;stage2 mbest;'); \
             quit" | octave-cli -p $(CODEC2)/octave

$(TRAIN)_lbg_b_res1.txt $(TRAIN)_b_lbg_res1.txt $(TRAIN)_b_lbg_mbest2.txt: $(TRAIN)_b.f32
	mbest=1 M=512 ../script/ratek_resampler.sh train_lbg train_120_b.f32
