

Modeling sample clock errors using sox:

sox -r 8000 -s -2 mod_dqpsk.raw -s -2 mod_dqpsk_8008hz.raw rate -h 8008

TODO

[X] Get file based mod and demod working again
[ ] timing wraps around
    + what is affect of bouncing back and forth over boundary?
    + could mean we are going back and forth a symbol
    + do we need to logic to lose or gain a frame?
    + think so, add or remove samples, or a whole frame
[ ] demod outputs ber (maybe after settling time)
[ ] try interfering sine wave
    + maybe swept
    + does modem fall over?
[ ] try non-flat channel, e.g. 3dB difference between hi and low tones
    + make sure all estimators keep working
[ ] test rx level sensitivity, i.e. 0 to 20dB attenuation
[ ] try to run from shell script
[ ] run a few tests
[ ] start coding in C and repeat tests
[ ] arb bit stream input to Octave mod and demod

Tests

[ ] fdmdv_demod('mod_dqpsk_8008hz.raw',1400*60);
[ ] fdmdv_demod('mod_dqpsk_7992hz.raw',1400*60);
[ ] mod_dqpsk_awgn_4dB_8008hz.raw
[ ] mod_dqpsk_good_4dB_8008hz.raw
[ ] mod_dqpsk_moderate_4dB_8008hz.raw
[ ] mod_dqpsk_moderate_4dB_7992hz.raw
