# ratek_resampler.mk

SHELL  := /bin/bash
CODEC2 := $(HOME)/codec2
TRAIN_FULL := ~/Downloads/train_120.spc
TRAIN      := train_120

ratek.png: $(TRAIN)_b_lbg_res1.txt $(TRAIN)_b_lbg_res2.txt
	echo "ratek_resampler_plot(\"ratek.png\", \
             \"$(TRAIN)_lbg_b_res1.txt\",'g-+;lbg1;', \
             \"$(TRAIN)_lbg_b_res2.txt\",'g-+;lbg2;'); \
             quit" | octave-cli -p $(CODEC2)/octave

$(TRAIN)_lbg_b_res1.txt $(TRAIN)_b_lbg_res1.txt: $(TRAIN)_b.f32
	M=512 ../script/ratek_resampler.sh train_lbg train_120_b.f32
