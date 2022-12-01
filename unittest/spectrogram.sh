#!/usr/bin/env bash
# spectrogram.sh
#
# Render a spectrogram from a wave file.

PATH=${PATH}:${HOME}/codec2/build_linux/src
CODEC2=${HOME}/codec2

fullfile=$1
filename=$(basename -- "$fullfile")
extension="${filename##*.}"
filename="${filename%.*}"

echo "pkg load signal; warning('off', 'all'); \
      s=load_raw('${fullfile}'); \
      plot_specgram(s, 8000, 500, 2500); print('${filename}.jpg', '-djpg'); \
      quit" | octave-cli -p ${CODEC2}/octave -qf > /dev/null
