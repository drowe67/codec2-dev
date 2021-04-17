#!/bin/bash
# ota_data.sh
#
# Automated Over The Air (OTA) data test for FreeDV OFDM HF data modems
#
# 1. Build codec2
# 2. Install kiwclient:
#    git clone git@github.com:jks-prv/kiwiclient.git
# 3. Install Hamlib cli tools
#
# TODO:
#  [ ] plot scatter diagram
#  [ ] SNR output
#  [ ] timestamps
#  [ ] option to listen to recording as it comes in

PATH=${PATH}:${HOME}/codec2/build_linux/src:${HOME}/kiwiclient

kiwi_url=""
port=8073
freq_kHz="7177"
tx_only=0
Nbursts=10
mode="datac0"

function print_help {
    echo
    echo "Automated Over The Air (OTA) data test for FreeDV OFDM HF data modems"
    echo
    echo "  usage ./ota_data.sh [-d] [-f freq_kHz] [-t] [-b Nbursts] [-p port] kiwi_url"
    echo
    echo "    -d       debug mode; trace script execution"
    echo "    -m mode  datac0|datac1|datac3"
    echo "    -t       Tx only, useful for manually observing SDRs which block multiple sessions from one IP"
    echo
    exit
}

POSITIONAL=()
while [[ $# -gt 0 ]]
do
key="$1"
case $key in
    -d)
        set -x	
        shift
    ;;
    -f)
        freq_kHz="$2"	
        shift
        shift
    ;;
    -m)
        mode="$2"	
        shift
        shift
    ;;
    -n)
        Nbursts="$2"	
        shift
        shift
    ;;
    -p)
        port="$2"	
        shift
        shift
    ;;
    -t)
        tx_only=1	
        shift
    ;;
    -h)
        print_help	
    ;;
    *)
    POSITIONAL+=("$1") # save it in an array for later
    shift
    ;;
esac
done
set -- "${POSITIONAL[@]}" # restore positional parameters

if [ $tx_only -eq 0 ]; then
    if [ $# -lt 1 ]; then
        print_help
    fi
    kiwi_url="$1"	
fi

# create test Tx file
freedv_data_raw_tx -q --framesperburst 1 --bursts ${Nbursts} --testframes ${Nbursts} ${mode} /dev/zero test_datac0.raw

usb_lsb=$(python3 -c "print('usb') if ${freq_kHz} >= 10000 else print('lsb')")

echo -n "waiting for KiwiSDR "
if [ $tx_only -eq 0 ]; then
    # start recording from remote kiwisdr
    kiwi_stdout=$(mktemp)
    kiwirecorder.py -s $kiwi_url -p ${port} -f $freq_kHz -m ${usb_lsb} -r 8000 --filename=rx --time-limit=60 >$kiwi_stdout &
    kiwi_pid=$!

    # wait for kiwi to start recording
    timeout_counter=0
    until grep -q -i 'Block: ' $kiwi_stdout
    do
        timeout_counter=$((timeout_counter+1))
        if [ $timeout_counter -eq 10 ]; then
            echo "can't connect to ${kiwi_url}"
            exit 1
        fi
        echo -n "."
        sleep 1
    done
    echo
fi

# transmit using local SSB radio
echo "Tx data signal"
freq_Hz=$((freq_kHz*1000))
usb_lsb_upper=$(echo ${usb_lsb} | awk '{print toupper($0)}')
echo "\\set_mode PKT${usb_lsb_upper} 0" | rigctl -m 361 -r /dev/ttyUSB0 > /dev/null
echo "\\set_freq ${freq_Hz}" | rigctl -m 361 -r /dev/ttyUSB0  > /dev/null
echo "\\set_ptt 1" | rigctl -m 361 -r /dev/ttyUSB0  > /dev/null
aplay --device="plughw:CARD=CODEC,DEV=0" -f S16_LE test_datac0.raw 2>/dev/null
echo "\\set_ptt 0" | rigctl -m 361 -r /dev/ttyUSB0  > /dev/null

if [ $tx_only -eq 0 ]; then
    sleep 2
    echo "Stopping KiwiSDR"
    kill ${kiwi_pid}
    wait ${kiwi_pid} 2>/dev/null

    echo "Process receiver sample"
    # generate spectrogram
    echo "pkg load signal; warning('off', 'all'); \
          s=load_raw('rx.wav'); \
          plot_specgram(s, 8000, 500, 2500); print('spec.png', '-dpng'); \
          quit" | octave-cli -p ../octave -qf
    # attempt to demodulate
    freedv_data_raw_rx -q --framesperburst 1 --testframes ${mode} -v --scatter scatter.txt rx.wav /dev/null
    # render scatter plot (if we get any frames)
    scatter_sz=$(ls -l scatter.txt | cut -f 5 -d' ')
    if [ $scatter_sz -ne 0 ]; then
        echo "pkg load signal; warning('off', 'all'); pl_scatter('scatter.txt'); quit"  | octave-cli -p ../octave -qf
    fi
fi
