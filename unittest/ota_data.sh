#!/bin/bash -x
# ota_data.sh
#
# Automated Over The Air (OTA) data test for FreeDV OFDM HF data modems
#
# 1. Build codec2
# 2. Install kiwclient:
#    git clone https://github.com/drowe67/freedv-gui.git

PATH=${PATH}:${HOME}/codec2/build_linux/src:${HOME}/kiwiclient
echo $PATH

# create test Tx file
freedv_data_raw_tx --framesperburst 1 --bursts 10 --testframes 10 datac0 /dev/zero test_datac0.raw

# start recording from remote kiwisdr
kiwirecorder.py -s sdrbris.proxy.kiwisdr.com -f 7177 -m lsb -r 8000 --filename=kiwirx --time-limit=12 &
kiwipid=$!
sleep 2

# transmit using local SSB radio
echo "\\set_ptt 1" | rigctl -m 361 -r /dev/ttyUSB0
aplay --device="plughw:CARD=CODEC,DEV=0" -f S16_LE test_datac0.raw
echo "\\set_ptt 0" | rigctl -m 361 -r /dev/ttyUSB0

wait ${kiwipid}

# attempt to demodulate
freedv_data_raw_rx --framesperburst 1 --testframes datac0 kiwirx.wav /dev/null
