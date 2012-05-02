

Modeling sample clock errors using sox:

sox -r 8000 -s -2 mod_dqpsk.raw -s -2 mod_dqpsk_8008hz.raw rate -h 8008

TODO

[ ] consider moving this file to root
    [ ] sep SVN repo, automake etc?
[ ] list each fdmdv.m script (ut, mod, demod) and C program/src and what it does
    [ ] example usage
[ ] repair fdmdv_ut and mod/demod after new statres exposed and var renames
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
    [ ] measure execution speed
[ ] document use of fdmdv_ut and fdmdv_demod + PathSim
[ ] more positibe form of sync reqd for DV frames?
    + like using track/acquire bit
[ ] more robust track/acquite state machine?
    + e.g. hang on thr fades?
[ ] block diagram
    [ ] maybe in ascii art
[ ] blog posts(s)
[ ] Codec 2 web page update
[ ] examples of various combinations of fdmdv demos
    $ ./fdmdv_get_test_bits test.c2 1400
    $ ./fdmdv_mod test.c2 test.raw
    $ play -r 8000 -s -2 test.raw
    
    Two seconds of test frame data modulated and sent out of sound device: 
       $ ./fdmdv_get_test_bits - 2800 | ./fdmdv_mod - - | play -t raw -r 8000 -s -2 -
 
    Count errors in two seconds of test frame data:
       $ ./fdmdv_get_test_bits - 2800 | ./fdmdv_put_test_bits -

    Ten sconds of modem simulation with testframes;

       $  ./fdmdv_get_test_bits - 14000 | ./fdmdv_mod - - | ./fdmdv_demod - - demod_dump.txt | ./fdmdv_put_test_bits -

[ ] PAPR idea
    + automatically tweak phases to reduce PAPR, e.g. slow variations in freq...
[ ] implement SNR and ppm measurements
[ ] add help to each octave script & C program

Tests

[ ] fdmdv_demod('mod_dqpsk_8008hz.raw',1400*60);
[ ] fdmdv_demod('mod_dqpsk_7992hz.raw',1400*60);
[ ] mod_dqpsk_awgn_4dB_8008hz.raw
[ ] mod_dqpsk_good_4dB_8008hz.raw
[ ] mod_dqpsk_moderate_4dB_8008hz.raw
[ ] mod_dqpsk_moderate_4dB_7992hz.raw
[ ] time ....
