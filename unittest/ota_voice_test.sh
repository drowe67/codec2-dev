#!/usr/bin/env bash
# ota_test.sh
#
# Automated Over The Air (OTA) voice test for FreeDV HF voice modes
#
# 1. Build codec2
# 2. Install kiwclient:
#    cd ~ && git clone git@github.com:jks-prv/kiwiclient.git
# 3. Install Hamlib cli tools

MY_PATH=`dirname $0`
BUILD_PATH=`echo $MY_PATH/../build_*/src`
PATH=${PATH}:${BUILD_PATH}:${HOME}/kiwiclient
CODEC2=${MY_PATH}/..

kiwi_url=""
port=8073
freq_kHz="7177"
tx_only=0
Nbursts=5
mode="700D"
model=361
gain=6
serialPort="/dev/ttyUSB0"
rxwavefile=0
soundDevice="plughw:CARD=CODEC,DEV=0"
txstats=0

function print_help {
    echo
    echo "Automated Over The Air (OTA) voice test for FreeDV HF voice modes"
    echo
    echo "  usage ./ota_voice_test.sh [options] SpeechFile [kiwi_url]"
    echo "  or:"
    echo "  usage ./ota_voice_test.sh -r rxWaveFile"
    echo
    echo "    -c dev    The sound device (in ALSA format on Linux, CoreAudio for macOS)"
    echo "    -d        debug mode; trace script execution"
    echo "    -g        SSB (analog) compressor gain"
    echo "    -m mode   700c|700d|700e"
    echo "    -o model  select radio model number ('rigctl -l' to list)"
    echo "    -r        Rx wave file mode: Rx process supplied rx wave file"
    echo "    -s port   The serial port (or hostname:port) to connect to for TX, default /dev/ttyUSB0"
    echo "    -t        Tx only, useful for manually observing SDRs"
    echo "    -x        Generate tx.wav file and exit"
    echo
    exit
}

# Approximation of Hilbert clipper type compressor.  Could do with some HF boost
function analog_compressor {
    input_file=$1
    output_file=$2
    gain=$3
    cat $input_file | ch - - 2>/dev/null | \
    ch - - --No -100 --clip 16384 --gain $gain 2>/dev/null | \
    # final line prints peak and CPAPR for SSB
    ch - - --clip 16384 |
    # manually adjusted to get similar peak levels for SSB and FreeDV
    sox -t .s16 -r 8000 -c 1 -v 0.85 - -t .s16 $output_file
}

function run_rigctl {
    command=$1
    model=$2
    echo $command | rigctl -m $model -r $serialPort > /dev/null
    if [ $? -ne 0 ]; then
        echo "Can't talk to Tx"
        exit 1
    fi
}

function clean_up {
    echo "killing KiwiSDR process"
    kill ${kiwi_pid}
    wait ${kiwi_pid} 2>/dev/null
    exit 1
}

function process_rx {
    echo "Process receiver sample"
    rx=$1
    # generate spectrogram
    echo "pkg load signal; warning('off', 'all'); \
          s=load_raw('${rx}'); \
          plot_specgram(s, 8000, 200, 3000); print('spec.jpg', '-djpg'); \
          quit" | octave-cli -p ${CODEC2}/octave -qf > /dev/null
    # attempt to decode
    freedv_rx ${mode} ${rx} - -v --highpassthroughgain 2>rx_stats.txt | sox -t .s16 -r 8000 -c 1 - rx_freedv.wav
    cat rx_stats.txt | tr -s ' ' | cut -f5 -d' ' | awk '$0==($0+0)' > sync.txt
    cat rx_stats.txt | tr -s ' ' | cut -f10 -d' ' | awk '$0==($0+0)' > snr.txt
    # time domain plot of output speech, SNR, and sync
    echo "pkg load signal; warning('off', 'all'); \
          s=load_raw('rx_freedv.wav'); snr=load('snr.txt'); sync=load('sync.txt'); \
          subplot(211); plot(s); subplot(212); x=1:length(sync); plotyy(x,snr,x,sync); \
          ylim([-5 15]); ylabel('SNR (dB)'); grid; \
          print('time_snr.jpg', '-djpg'); \
          printf('Nsync: %3d\n', sum(sync)); \
          snr_valid = snr(find(snr != -5.0)); \
          if length(snr_valid) printf('SNRav: %5.2f\n', mean(snr_valid)); else printf('SNRav: %5.2f\n', -5); end;
          quit" | octave-cli -p ${CODEC2}/octave -qf
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
    -g)
        gain="$2"	
        shift
        shift
    ;;
    -o)
        model="$2"	
        shift
        shift
    ;;
    -m)
        mode="$2"	
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
    -r)
        rxwavefile=1	
        shift
    ;;
    -x)
        txstats=1	
        shift
    ;;
    -c)
        soundDevice="$2"
        shift
        shift
    ;;
    -s)
        serialPort="$2"
        shift
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

if [ $rxwavefile -eq 1 ]; then
    process_rx $1
    exit 0
fi

speechfile="$1"
if [ ! -f $speechfile ]; then
    echo "Can't find ${speechfile}!"
    exit 1
fi

if [ $tx_only -eq 0 ]; then
    if [ $# -lt 1 ]; then
        print_help
    fi
    kiwi_url="$2"
    echo $kiwi_url
fi

# create Tx file ------------------------
echo $mode

# create compressed analog
speech_comp=$(mktemp)
speech_freedv=$(mktemp)
analog_compressor $speechfile $speech_comp $gain

# create modulated FreeDV, with compressor enabled
freedv_tx $mode $speechfile $speech_freedv --clip 1 
cat $speech_comp $speech_freedv > tx.raw
sox -t .s16 -r 8000 -c 1 tx.raw tx.wav

if [ $txstats -eq 1 ]; then
    # ch just used to monitor observe peak and RMS level
    ch $speech_freedv /dev/null
    # time domain plot of tx signal
    echo "pkg load signal; warning('off', 'all'); \
          s=load_raw('tx.raw'); plot(s); \
          print('tx.jpg', '-djpg'); \
          quit" | octave-cli -p ${CODEC2}/octave -qf > /dev/null
    exit 0
fi

# kick off KiwiSDR ----------------------------

usb_lsb=$(python3 -c "print('usb') if ${freq_kHz} >= 10000 else print('lsb')")
if [ $tx_only -eq 0 ]; then
    # clean up any kiwiSDR processes if we get a ctrl-C
    trap clean_up SIGHUP SIGINT SIGTERM

    echo -n "waiting for KiwiSDR "
    # start recording from remote kiwisdr
    kiwi_stdout=$(mktemp)
    kiwirecorder.py -s $kiwi_url -p ${port} -f $freq_kHz -m ${usb_lsb} -r 8000 --filename=rx --time-limit=300 >$kiwi_stdout &
    kiwi_pid=$!

    # wait for kiwi to start recording
    timeout_counter=0
    until grep -q -i 'Block: ' $kiwi_stdout
    do
        timeout_counter=$((timeout_counter+1))
        if [ $timeout_counter -eq 10 ]; then
            echo "can't connect to ${kiwi_url}"
            kill ${kiwi_pid}
            wait ${kiwi_pid} 2>/dev/null
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
run_rigctl "\\set_mode PKT${usb_lsb_upper} 0" $model
run_rigctl "\\set_freq ${freq_Hz}" $model
run_rigctl "\\set_ptt 1" $model
if [ `uname` == "Darwin" ]; then
    AUDIODEV="${soundDevice}" play -t raw -b 16 -c 1 -r 8000 -e signed-integer --endian little tx.raw 
else
    aplay --device="${soundDevice}" -f S16_LE tx.raw 2>/dev/null
fi
run_rigctl "\\set_ptt 0" $model

if [ $tx_only -eq 0 ]; then
    sleep 2
    echo "Stopping KiwiSDR"
    kill ${kiwi_pid}
    wait ${kiwi_pid} 2>/dev/null

    process_rx rx.wav
fi

