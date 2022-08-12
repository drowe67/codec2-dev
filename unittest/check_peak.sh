#!/usr/bin/env bash
#
# Check peak level of each FreeDV waveform is about the same to present
# consistent drive to transmitters.
#
# For manual run outside of ctest:
#  cd codec/build_linux
#  PATH=${PATH}:${HOME}/codec2/build_linux/src
#  ./unittest/check_peak.sh
#  OR:
#  

voice_test() {
    mode=$1
    echo -n "$mode "
    f=$(mktemp)
    freedv_tx $mode ../raw/ve9qrp_10s.raw $f --clip 1
    octave_cmd="cd ../octave; 
                t=load_raw('${f}'); 
                mx=max(t); printf('%d ',max(t)); 
                if (mx > 16000) && (mx < 17000) printf('PASS\n') else printf('FAIL\n') end"
    octave-cli -qf --eval "$octave_cmd"
}

data_test() {
    mode=$1
    echo -n "$mode "
    f=$(mktemp)
    freedv_data_raw_tx --framesperburst 2 --bursts 3 --testframes 6 $mode /dev/zero $f 2>/dev/null
    octave_cmd="cd ../octave; 
                t=load_raw('${f}'); 
                mx=max(t); printf('%d ',max(t)); 
                if (mx > 16000) && (mx < 17000) printf('PASS\n') else printf('FAIL\n') end"
    octave-cli -qf --eval "$octave_cmd"
}

if [ "$1" == "LPCNet" ]; then
    # these don't get run unless we build with LPCNet
    voice_test "2020"
    voice_test "2020B"
    voice_test "2020C"
  else
    voice_test "1600"
    voice_test "700C"
    voice_test "700D"
    voice_test "700E"
    voice_test "800XA"
    voice_test "2400A"
    voice_test "2400B"
    data_test "datac0"
    data_test "datac1"
    data_test "datac3"
fi

exit 0

