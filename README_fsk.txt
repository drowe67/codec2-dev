
README_fsk.txt
David Rowe
Created Jan 2016

A FSK modem with a non-coherent demodulator.  Performance is within a
fraction of a dB of ideal.  The demodulator can estimate and track
frequency offsets.

Credits
-------

The Octave version of the modem was developed by David Rowe.  Brady
O'Brien ported the modem to C, and wrote the C/Octave tests.

Quickstart
----------

Built as part of codec2-dev, see README for build instructions.

1. Using a random bit stream input generate 100 bit/s FSK and play out
   speaker:
    
    $ cd build_linux/src
    $ /fsk_get_test_bits - 10 | ./fsk_mod 8000 100 1200 2400 - - | play -t raw -r 8000 -s -2 

2. Check the bit error rate of 100 frames of 1200 bit/s FSK:

    $ ./fsk_get_test_bits - 100 | ./fsk_mod 9600 1200 1200 2400 - - | ./fsk_demod 9600 1200 - - | ./fsk_put_test_bits -
      FSK BER 0.000000, bits tested 9499, bit errors 0

3. Automatically check C implementation against Octave using a variety of test conditions and 
   channel impairments.

    $ cd octave
    $ octave
    octave:1> tfsk

