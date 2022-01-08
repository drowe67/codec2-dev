#!/usr/bin/env bash
#
# Summarise tests to date into one directory to allow easy browsing

dir=voice_summary
mkdir -p ${dir}
time_snr_files=$(find . -name time_snr.jpg | sort)
p=$(pwd)
serial=0
for f in $time_snr_files
do
    d=$(echo $f | sed -r 's/\.\/(.*)\/time_snr.jpg/\1/')
    sdr_url=$(head ${d}/log.txt -n 1)
    sdr="unk"
    case $sdr_url in
        "kiwisdr.areg.org.au")
            sdr="areg"
            ;;
        "sdr-amradioantennas.com")
            sdr="am"
            ;;
        "vk6qs.proxy.kiwisdr.com")
            sdr="vk6qs"
            ;;
        "sdr.ironstonerange.com")
            sdr="iron"
            ;;
        "kk6pr.ddns.net")
            sdr="kk6pr"
            ;;
        "kiwisdr.owdjim.gen.nz")
            sdr="marahau"
            ;;
        "kiwisdrzl1kfm.ddns.net")
            sdr="zl1kfm"
            ;;
        *)
            echo "Unknown Kiwi SDR"
            ;;
    esac
    mode=$(head ${d}/log.txt -n 2 | tail -n 1)
    serial_str=$(printf "%04d" $serial)
    echo $serial_str $d $sdr $mode
    cp ${d}/spec.jpg ${dir}/${serial_str}_${d}_${sdr}_${mode}_spec.jpg
    cp ${d}/time_snr.jpg ${dir}/${serial_str}_${d}_${sdr}_${mode}_time_snr.jpg
    cp ${d}/rx.wav ${dir}/${serial_str}_${d}_${sdr}_${mode}_rx.wav
    cp ${d}/rx_freedv.wav ${dir}/${serial_str}_${d}_${sdr}_${mode}_rx_freedv.wav
    serial=$((serial + 1))
done   
