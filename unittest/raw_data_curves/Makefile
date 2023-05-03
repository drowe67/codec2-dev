# Makefile
# Dec 2022
#
# Automates PER/BER curve generation for raw data mode:
#
# 1. Compare "ch" noise injection/SNR measurement against reference Octave Tx
# 2. Compare C Tx against Octave Tx (with and without compression)
# 3. Plot curves for SNR estimated from C Rx against actual SNR
# 4. Plot AWGN/PER C Tx curves for end user documentation

SHELL  := /bin/bash
CODEC2 := $(HOME)/codec2

all: test \
     octave_ch_noise_awgn.png octave_c_tx_awgn.png octave_c_tx_comp_awgn.png \
     octave_ch_noise_mpp.png octave_c_tx_mpp.png octave_c_tx_comp_mpp.png \
     snrest_snr_ctx.png snrest_snr_ctxc.png \
     c_tx_comp.png c_tx_comp_thruput.png

clean:
	rm -f *.txt *.png *.raw

# run this first, traps common CML setup error
test:
	source snr_curves.sh; test_ldpc

# subset of files generated, but enough to set up Makefile dependencies
snr_oct =  snr_oct_datac0_awgn.txt  snr_oct_datac1_awgn.txt  snr_oct_datac3_awgn.txt
snr_ch  =  snr_ch_datac0_awgn.txt   snr_ch_datac1_awgn.txt   snr_ch_datac3_awgn.txt
snr_ctx =  snr_ctx_datac0_awgn.txt  snr_ctx_datac1_awgn.txt  snr_ctx_datac3_awgn.txt
snr_ctxc = snr_ctxc_datac0_awgn.txt snr_ctxc_datac3_awgn.txt

snr_oct_mpp  = snr_oct_datac0_mpp.txt  snr_oct_datac1_mpp.txt  snr_oct_datac3_mpp.txt
snr_ch_mpp   = snr_ch_datac0_mpp.txt   snr_ch_datac1_mpp.txt   snr_ch_datac3_mpp.txt
snr_ctx_mpp  = snr_ctx_datac0_mpp.txt  snr_ctx_datac1_mpp.txt  snr_ctx_datac3_mpp.txt
snr_ctxc_mpp = snr_ctxc_datac0_mpp.txt snr_ctxc_datac3_mpp.txt

$(snr_oct):
	source snr_curves.sh; generate_octave_tx_data datac0 awgn
	source snr_curves.sh; generate_octave_tx_data datac1 awgn
	source snr_curves.sh; generate_octave_tx_data datac3 awgn

$(snr_oct_mpp):
	source snr_curves.sh; generate_octave_tx_data datac0 mpp
	source snr_curves.sh; generate_octave_tx_data datac1 mpp
	source snr_curves.sh; generate_octave_tx_data datac3 mpp

$(snr_ch):
	source snr_curves.sh; generate_ch_data datac0 awgn
	source snr_curves.sh; generate_ch_data datac1 awgn
	source snr_curves.sh; generate_ch_data datac3 awgn

$(snr_ch_mpp):
	source snr_curves.sh; generate_ch_data datac0 mpp
	source snr_curves.sh; generate_ch_data datac1 mpp
	source snr_curves.sh; generate_ch_data datac3 mpp

# C without compression

$(snr_ctx):
	source snr_curves.sh; generate_snrest_v_snr_data datac0 awgn
	source snr_curves.sh; generate_snrest_v_snr_data datac1 awgn
	source snr_curves.sh; generate_snrest_v_snr_data datac3 awgn
	source snr_curves.sh; generate_snrest_v_snr_data datac4 awgn
	source snr_curves.sh; generate_snrest_v_snr_data datac13 awgn

$(snr_ctx_mpp):
	source snr_curves.sh; generate_snrest_v_snr_data datac0 mpp
	source snr_curves.sh; generate_snrest_v_snr_data datac1 mpp
	source snr_curves.sh; generate_snrest_v_snr_data datac3 mpp
	source snr_curves.sh; generate_snrest_v_snr_data datac4 mpp
	source snr_curves.sh; generate_snrest_v_snr_data datac13 mpp

# C with compression

$(snr_ctxc):
	source snr_curves.sh; generate_snrest_v_snr_data datac0 awgn 1
	source snr_curves.sh; generate_snrest_v_snr_data datac1 awgn 1
	source snr_curves.sh; generate_snrest_v_snr_data datac3 awgn 1
	source snr_curves.sh; generate_snrest_v_snr_data datac4 awgn 1
	source snr_curves.sh; generate_snrest_v_snr_data datac13 awgn 1

$(snr_ctxc_mpp):
	source snr_curves.sh; generate_snrest_v_snr_data datac0 mpp 1
	source snr_curves.sh; generate_snrest_v_snr_data datac1 mpp 1
	source snr_curves.sh; generate_snrest_v_snr_data datac3 mpp 1
	source snr_curves.sh; generate_snrest_v_snr_data datac4 mpp 1
	source snr_curves.sh; generate_snrest_v_snr_data datac13 mpp 1

# Octave and C curves should be on top of each other, indicating Octave
# and ch noise injection/SNR measurement are equivalent
octave_ch_noise_awgn.png: $(snr_oct) $(snr_ch)
	echo "snr_curves_plot; octave_ch_noise_print('awgn'); quit" | \
        octave-cli -p $(CODEC2)/octave
octave_ch_noise_mpp.png: $(snr_oct_mpp) $(snr_ch_mpp)
	echo "snr_curves_plot; octave_ch_noise_print('mpp'); quit" | \
        octave-cli -p $(CODEC2)/octave

# Octave Tx and C Tx curves should be on top of each other
octave_c_tx_awgn.png: $(snr_oct) $(snr_ctx)
	echo "snr_curves_plot; octave_c_tx_print('awgn'); quit" | \
        octave-cli -p $(CODEC2)/octave
octave_c_tx_mpp.png: $(snr_oct_mpp) $(snr_ctx_mpp)
	echo "snr_curves_plot; octave_c_tx_print('mpp'); quit" | \
        octave-cli -p $(CODEC2)/octave

# Octave Tx and C Tx (compressed) curves should be close, but C may be 1dB
# poorer
octave_c_tx_comp_awgn.png: $(snr_oc) $(snr_ctxc)
	echo "snr_curves_plot; octave_c_tx_comp_print('awgn'); quit" | \
        octave-cli -p $(CODEC2)/octave
octave_c_tx_comp_mpp.png: $(snr_oct_mpp) $(snr_ctxc_mpp)
	echo "snr_curves_plot; octave_c_tx_comp_print('mpp'); quit" | \
        octave-cli -p $(CODEC2)/octave

# combined AWGN and MPP from C Tx (compressed) - what end users would run
c_tx_comp.png: $(snr_ctxc) $(snr_ctxc_mpp)
	echo "snr_curves_plot; c_tx_comp_print; quit" | \
        octave-cli -p $(CODEC2)/octave

# Curves of SNR estimates from C Rx compared to actual SNR, useful for "gear shifting"
snrest_snr_ctx.png: $(snr_ctx)
	echo "snr_curves_plot; snrest_snr_print('ctx', 'awgn'); quit" | \
        octave-cli -p $(CODEC2)/octave
snrest_snr_ctxc.png: $(snr_ctxc)
	echo "snr_curves_plot; snrest_snr_print('ctxc', 'awgn'); quit" | \
        octave-cli -p $(CODEC2)/octave

# Throughput of payload data in bits/s of modes against SNR
c_tx_comp_thruput.png: $(snr_ctxc) $(snr_ctxc_mpp)
	echo "snr_curves_plot; c_tx_comp_thruput_print; quit" | \
        octave-cli -p $(CODEC2)/octave
