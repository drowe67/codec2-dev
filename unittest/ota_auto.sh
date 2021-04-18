#!/bin/bash -x
# ota_auto.sh
#
# Run a single automated test and log results

timestamp=$(date +"%F-%T")
mkdir -p $timestamp
start_dir=$(pwd)
cd $timestamp
../ota_data.sh "$@" >> log.txt 2>&1
cd $start_dir
kiwi_sdr=$(head -n 1 ${timestamp}/log.txt)
result=$(tail -n 1 ${timestamp}/log.txt)
echo $timestamp $kiwi_sdr $result >> log.txt
