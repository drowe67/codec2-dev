
README_fsk.txt
David Rowe
Created Jan 2016

A FSK modem with a non-coherent demodulator.  Performance is within a
fraction of a dB of ideal.  The demodulator can estimate the tone
frequencies and track drift.

Credits
-------

The Octave version of the modem was developed by David Rowe.  Brady
O'Brien ported the modem to C, and wrote the C/Octave tests.

Quickstart
----------

Built as part of codec2-dev, see README for build instructions.

1. Using 10 test frames as a source of bits, generate 2FSK using a
   8000 Hz sample rate, at 100 symbols/s (== 100 bits/ for 2FSK), with
   two frequencies of 1200 and 2400 Hz, and stream via stdout to
   "play" out speaker:
    
    $ cd build_linux/src
    $ ./fsk_get_test_bits - 10 | ./fsk_mod 2 8000 100 1200 2400 - - | play -t raw -r 8000 -s -2 -

    Notes:
      i) TODO explain structure of test frame and how many bits are in it.

2. Measure the bit error rate of 100 frames of 100 bit/s 2FSK:

    $ ./fsk_get_test_bits - 100 | ./fsk_mod 2 8000 100 1200 2400 - - | ./fsk_demod 2 0 8000 100 - - | ./fsk_put_test_bits -
       FSK BER 0.003108, bits tested 39899, bit errors 124

    Notes:
      i) The "P" parameter is don't care when mode is 2 or 4, so we set it to 0
      ii) TODO: waht are the source of errors?  Can we remove them with no noise?

3. Measure the bit error rate of 100 frames at 1200 bits/s, samples at 9600 Hz:

    $ ./fsk_get_test_bits - 100 | ./fsk_mod 2 9600 1200 1200 2400 - - | ./fsk_demod 2X 8 9600 1200 - - | ./fsk_put_test_bits -
      FSK BER 0.000000, bits tested 9499, bit errors 0

    Notes: i) The 2X option is used in the demod, which narrows the
              window for parameter estimation window.  In 2X/4X mode
              the P parameter is used so we set it to 8 (TODO: explain why)
           ii) TODO explain Restrictions to bit rate ie integer P, so we changed Fs to suit?

4.  (TODO High bit rate example like 115k project Horus)


5. Automatically check C implementation against Octave using a variety of test conditions and 
   channel impairments.

    $ cd octave
    $ octave
    octave:1> tfsk

