
README_fsk.txt
David Rowe
Created Jan 2016

A FSK modem with a non-coherent demodulator.  Performance is within a
fraction of a dB of ideal.  The demodulator automagically estimates
the tone frequencies and tracks frequency drift.

Credits
-------

The Octave version of the modem was developed by David Rowe.  Brady
O'Brien ported the modem to C, and wrote the C/Octave tests.

Quickstart
----------

Built as part of codec2-dev, see README for build instructions.

1. Using 10 test frames as a source of bits, generate 2FSK using a
   8000 Hz sample rate, at 100 symbols/s (== 100 bit/s for 2FSK), with
   two frequencies of 1200 and 2400 Hz, and stream via stdout to
   "play" out speaker:
    
    $ cd build_linux/src
    $ ./fsk_get_test_bits - 10 | ./fsk_mod 2 8000 100 1200 1200 - - | play -t raw -r 8000 -s -2 -

    Each test frame has 400 bits, and the "fsk_get_test_bits"
    generates 1+10=11 test frames in the example above - an extra
    first frame is always generated to give the demod time to sync.

    The low tone frequency is 1200Hz, and the upper tone 1200 + 1200 = 2400Hz.

2. Measure the bit error rate of 100 frames of 100 bit/s 2FSK:

    $ ./fsk_get_test_bits - 100 | ./fsk_mod 2 8000 100 1200 100 - - | ./fsk_demod 2 0 8000 100 - - | ./fsk_put_test_bits -
    FSK BER 0.000000, bits tested 39899, bit errors 0

    In this example the two tones are at 1200 and 1200+100 = 1300Hz.
    A shift of 100Hz is the minimum possible for an incoherent FSK
    demodulator that is running at 100 symbols/s.

    Note the "P" parameter (explained below) is don't care when mode
    is 2 or 4, so we set it to 0.

3. Measure the bit error rate of 100 frames at 1200 bits/s, using a
   sample rate of 9600 Hz:

    $ ./fsk_get_test_bits - 100 | ./fsk_mod 2 9600 1200 1200 1200 - - | ./fsk_demod 2X 8 9600 1200 - - | ./fsk_put_test_bits -
    FSK BER 0.000000, bits tested 39967, bit errors 0

    The 2X option is a "high speed" mode, which narrows the parameter
    estimation window to make demod synchronise quickly in low-latency
    applications like PTT digital voice.  For low speed telemetry
    applications we use low speed mode (a wider window) as processing
    latency isn't important and we want the best possible estimate.
   
    In 2X/4X mode the P parameter _is_ used, in this example we set it
    to 8.  This means the timing estimator uses a window that is 8
    symbols wide.  There are some restrictions on P, Fs/(Rs*P) needs
    to be an integer.  So for Fs=9600Hz, Rs=1200Hz Fs/Rs=8, so P=8
    would be OK but P=10 not OK.  An assert will fire if you get it
    wrong.

4.  (TODO High bit rate example like 115k project Horus)

5. Automatically check C implementation against Octave using a variety of test conditions and 
   channel impairments.

    $ cd octave
    $ octave
    octave:1> tfsk

