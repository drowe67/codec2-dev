

Modeling sample clock errors using sox:

sox -r 8000 -s -2 mod_dqpsk.raw -s -2 mod_dqpsk_8008hz.raw rate -h 8008

TODO

[X] Get file based mod and demod working again
[ ] try interfering sine wave
    + maybe swept
    + does modem fall over?
[ ] try non-flat channel, e.g. 3dB difference between hi and low tones
    + make sure all estimators keep working
[ ] test rx level sensitivity, i.e. 0 to 20dB attenuation
[ ] try to run from shell script
[ ] arb bit stream input to Octave mod and demod
[ ] make fine freq indep of amplitude
    + use angle rather than imag corrd
[ ] C port and UT framework
    [ ] document how to use
    [ ] demo modem C test program
    [ ] freq corr in func, state vars in struct
    [ ] fine freq est in func, statevars
    [ ] demod in func with all vars
    [ ] mod in func with all vars
    [ ] check with ch impairments    
    [ ] test with freq offsets
[ ] document use of fdmdv_ut and fdmdv_demod + PathSim
[ ] block diagram
    [ ] maybe in ascii art
[ ] blog posts(s)
[ ] Codec 2 web page update

Tests

[ ] fdmdv_demod('mod_dqpsk_8008hz.raw',1400*60);
[ ] fdmdv_demod('mod_dqpsk_7992hz.raw',1400*60);
[ ] mod_dqpsk_awgn_4dB_8008hz.raw
[ ] mod_dqpsk_good_4dB_8008hz.raw
[ ] mod_dqpsk_moderate_4dB_8008hz.raw
[ ] mod_dqpsk_moderate_4dB_7992hz.raw
