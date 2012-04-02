

Modeling sample clock errors using sox:

sox -r 8000 -s -2 mod_dqpsk.raw -s -2 mod_dqpsk_8008hz.raw rate -h 8008

TODO

[ ] Get file based mod and demod working again
[ ] write file based ch simulator
[ ] demod outputs ber (maybe after settling time)
[ ] try to run from shell script
[ ] run a few tests
[ ] start coding in C and repeat tests
[ ] arb bit stream input to Octave mod and demod

