# Makefile for raw_data_curves
# Dec 2022

SHELL  := /bin/bash
CODEC2 := $(HOME)/codec2

all: octave_ch_noise.png octave_c_tx.png octave_c_tx_comp.png snrest_snr.png

clean:
	rm -f *.txt *.png *.raw

snr_oct = snr_oct_datac0.txt snr_oct_datac1.txt snr_oct_datac3.txt
snr_ch  = snr_ch_datac0.txt snr_ch_datac1.txt snr_ch_datac3.txt
snr_ctx = snr_ctx_datac0.txt snr_ctx_datac1.txt snr_ctx_datac3.txt
snr_ctxc = snr_ctxc_datac0.txt snr_ctxc_datac3.txt

$(snr_oct):
	source snr_curves.sh; generate_octave_tx_data 'datac0'
	source snr_curves.sh; generate_octave_tx_data 'datac1'
	source snr_curves.sh; generate_octave_tx_data 'datac3'

$(snr_ch):
	source snr_curves.sh; generate_ch_data 'datac0'
	source snr_curves.sh; generate_ch_data 'datac1'
	source snr_curves.sh; generate_ch_data 'datac3'

$(snr_ctx):
	source snr_curves.sh; generate_snrest_v_snr_data 'datac0'
	source snr_curves.sh; generate_snrest_v_snr_data 'datac1'
	source snr_curves.sh; generate_snrest_v_snr_data 'datac3'

$(snr_ctxc):
	source snr_curves.sh; generate_snrest_v_snr_data 'datac0' 1
	source snr_curves.sh; generate_snrest_v_snr_data 'datac3' 1

octave_ch_noise.png: $(snr_oct) $(snr_ch)
	echo "snr_curves_plot; octave_ch_noise_print; quit" | \
        DISPLAY="" octave-cli -p $(CODEC2)/octave

octave_c_tx.png: $(snr_oct) $(snr_ctx)
	echo "snr_curves_plot; octave_c_tx_print; quit" | \
        DISPLAY="" octave-cli -p $(CODEC2)/octave

octave_c_tx_comp.png: $(snr_oct) $(snr_ctxc)
	echo "snr_curves_plot; octave_c_tx_comp_print; quit" | \
        DISPLAY="" octave-cli -p $(CODEC2)/octave

snrest_snr.png: $(snr_ctx)
	echo "snr_curves_plot; snrest_snr_print; quit" | \
        DISPLAY="" octave-cli -p $(CODEC2)/octave
