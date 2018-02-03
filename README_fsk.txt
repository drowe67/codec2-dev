
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

1. Using 1000 test bits, generate 2FSK using a 8000 Hz sample rate, at
   100 symbols/s (== 100 bit/s for 2FSK), with two frequencies of 1200
   and 2400 Hz, and stream via stdout to "play" out speaker:
    
    $ cd build_linux/src
    $ ./fsk_get_test_bits - 1000 | ./fsk_mod 2 8000 100 1200 1200 - - | play -t raw -r 8000 -s -2 -

    The low tone frequency is 1200Hz, and the upper tone 1200 + 1200 = 2400Hz.

2. Measure the bit error rate over 10,000 bits of 100 bit/s 2FSK:

    $ ./fsk_get_test_bits - 10000 | ./fsk_mod 2 8000 100 1200 100 - - | ./fsk_demod -l 2 8000 100 - - | ./fsk_put_test_bits -
    FSK BER 0.000000, bits tested 9800, bit errors 0

    In this example the two tones are at 1200 and 1200+100 = 1300Hz.
    A shift of 100Hz is the minimum possible for an incoherent FSK
    demodulator that is running at 100 symbols/s. Note that -l or --lbr
    initalizes the demodulator in 'low bit rate' mode. In low bit rate
    mode, incoming samples are processed in 1 second chunks, giving a 
    very wide window for frequency and timing estimaton. Low speed mode
    is well suited to applications that can tolerate long latency, such
    as balloon telemetry.

    The demod and test frame logic takes a few frames to sync up, so
    although we sent 10,000 bits, only 9800 are received.  However
    there were no errors in those received bits.

3. Measure the bit error rate of 5000 bits at 1200 bits/s, using a
   sample rate of 9600 Hz:

    $ ./fsk_get_test_bits - 5000 | ./fsk_mod 2 9600 1200 1200 1200 - - | ./fsk_demod -p 8 2 9600 1200 - - | ./fsk_put_test_bits -
    FSK BER 0.000000, bits tested 4900, bit errors 0

    In this example, the -l and --lbr options are left out setting the
    modem up in "high speed" mode. In this mode, the demodulator
    operates on blocks of 24 symbols at a time. This allows the modem
    to operate with a lower latency. High speed mode is well suited to
    applications without much tolerance for latency or with limited
    processing power and memory, such as PTT digital voice.
    
    In this example, the -p (or --conv) option is used. This specifies
    the downconverted symbol size. The symbol period (Ts) must be
    divisible by the supplied P parameter. In this case, Fs is 9600
    and Rs is 1200, so Ts is Fs/Rs, which is 8. In this case, 8 is the
    maximum P value allowed, though 4 and 2 would also work. If P is
    not divisible by Ts, a failed assert will fire and fsk_demod will
    exit. If -p is not supplied, it will default to match Ts. In
    general, lower values of P result in less memory use and CPU, but
    may result in worse modem preformance.  P should be left alone
    unless CPU and memory usage needs to be lowered for the
    application.

    (TODO, make this easier to understand, perhaps with figure)

4.  (TODO High bit rate example like 115k project Horus)

5. Automatically check C implementation against Octave using a variety of test conditions and 
   channel impairments.

    $ cd octave
    $ octave
    octave:1> tfsk

