# FreeDV README

FreeDV is an open source digital voice protocol that integrates modems, codecs, and FEC. 

On transmit, FreeDV converts speech to a modem signal you can send over a radio channel.  On receive, FreeDV takes off air modem signals and converts them to speech samples.

FreeDV is available as a GUI application, an open source library (FreeDV API), and in hardware (the SM1000 FreeDV adaptor).  FreeDV is part of the Codec 2 project.

![FreeDV mode knob](http://www.rowetel.com/images/codec2/mode_dv.jpg)

## Reading Further

1. [FreeDV web site](http://freedv.org)
1. [FreeDV GUI User Manual](https://github.com/drowe67/freedv-gui/blob/master/USER_MANUAL.md)
1. [Codec 2](http://rowetel.com/codec2.html)
1. FreeDV can also be used for data [README_data](https://github.com/drowe67/codec2/blob/master/README_data.md)
1. [FreeDv 700C blog post](http://www.rowetel.com/wordpress/?p=5456)
1. [FreeDV 700D Released blog post](http://www.rowetel.com/wordpress/?p=6103)
1. [FreeDV 2020 blog post](http://www.rowetel.com/wordpress/?p=6747)
1. [FreeDV 2400A blog post](http://www.rowetel.com/?p=5119)
1. [FreeDV 2400A & 2400B](http://www.rowetel.com/?p=5219)
1. Technical information on various modem waveforms in the [Waveform spreadsheet](https://github.com/drowe67/codec2/blob/master/doc/modem_codec_frame_design.ods)

## FreeDV HF Modes

These are designed for use with a HF SSB radio.

| Mode | Date | Codec2 Mode | Modem | Raw bits/s | FEC | Text | SNR min | Multipath |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1600 | 2012 | 1300 | DQPSK + pilot | 1600 | Golay (23,12) | Y | 2 | poor |
| 700C | 2017 | 700C | coherent QPSK + diversity | 1400 | - | - | 2 | good |
| 700D | 2018 | 700C | coherent OFDM/QPSK | 1900 | LDPC (224,112)  | Y | -2 | fair |
| 2020 | 2019 | LPCNet | coherent OFDM/QPSK | 3000 | LDPC (504,396)  | Y | 2 | poor |

TODO RF bandwidth

Notes:

1. All of the modems use multiple parallel carriers running at a low symbol rate of aorund 50 Hz.  This helps combat the effects of multipath channels.

1. Coherent demodulation gives much better performance than differential, at the cost of some additional compelxity.  Pilot symbols are transmitted regularaly to allow the demod to estimate the channels phasemodems perform better that differential.  Additional complexity.

1. The 1600 and 700C waveforms use parallel tone modems, later modes use OFDM.  OFDM gives tighter carrier packing which allows higher bit rates, but tends to suffer more from frequency offsets and delay spread.

1. Raw bits/s is the number of payload bits/s carried over the channel by the modem.  This consists of codec frames, FEC parity bits, unprotected text, and synchronisation informations such as pilot and unique word bits.  The estimates are opento interpretation for the OFDM waveforms due to pilot symbol and cyclic prefix considerations (see spreadsheet).

1. At medium to high SNRs, FreeDV 700C performs well (better than 700D) on fast fading multipath channels with large delay spread, due it's parallel tone design, and high pilot symbol rate.  It employs transmit diversity which delivers BER performance similar to modes using FEC.  FreeDV 700C also has a short frame (40ms), so syncs fast with low latency.  This is useful on marginal channels that move between unsuable and barely usable.

1. FreeDV 700D uses an OFDM modem and was optimised for low SNR channels, with strong FEC but a low pilot symbol rate and modest (2ms) cyclic prefix which means it's performance degrades on multipath channels with fast (> 1Hz) fading.  The use of strong FEC makes this mode quite robust to other channel impairments, such as static crashes, urban HF noise, and in-band interference.

1. The "Text" is a side channel for low bit rate text such as your location and callsign.  It is generally unprotected, and encoded with varicode.

1. FEC was added fairly recently to FreeDV modes.  The voice codecs we use work OK at bit error rates of a few %, and packet error rates of 10%. Raw bit error rates on multipath channels often exceed 10%.  Thus to be useful we require a FEC code that works at over 10% raw BER, and has 1% output (coded) bit error rate.  This is an unusual requirement, most FEC codes are designed for data which is intolerant of any bit errors.  However powerful LDPC codes come close, and can also "clean up" other channel imparements such as static and interference.  The use of OFDM means we now have "room" for the extra bits required for FEC, so there is no cost in adding it.

## FreeDV VHF Modes

The modes used FSK or FM.

| Mode | Date | Codec2 Mode | Modem | Raw bits/s | FEC | Text |
| --- | --- | --- | --- | --- | --- | --- |
| 2400A | 2016 | 1300 | 4FSK | 2400 | Golay (23,12) | Y |
| 2400B | 2016 | 1300 | baseband/analog FM | 2400 | Golay (23,12) | Y |
| 800XA | 2017 | 700C |  4FSK | 800 | - | N |

## FreeDV Command line Software

```freedv_tx``` and ```freedv_rx``` are command line implementations of the FreeDV protocol, which combines Codec 2, modems, and Forward Error Correction (FEC).
  
## FreeDV API

See ```freedv_api.h``` and ```freedv_api.c```, and the demo programs ```freedv_tx``` & ```freedv_rx``.  Quickstart demo using FreeDV 1600:
```
$ ./freedv_tx 1600 ../../raw/hts1.raw - | ./freedv_rx 1600 - - | play -t raw -r 8000 -s -2 -q -
$ cat freedv_rx_log.txt
```
## FreeDV 2400A and 2400B modes

FreeDV 2400A and FreeDV 2400B are modes designed for VHF radio. FreeDV 2400A is designed for SDR radios (it has a 5 kHz RF bandwidth), however FreeDV 2400B is designed to pass through commodity FM radios.

Demos of FreeDV 2400A and 2400B:
```
$ ./freedv_tx 2400A ../../raw/ve9qrp_10s.raw - | ./freedv_rx 2400A - - | play -t raw -r 8000 -s -2 -
$ ./freedv_tx 2400B ../../raw/ve9qrp_10s.raw - | ./freedv_rx 2400B - - | play -t raw -r 8000 -s -2 -
```
Note for FreeDV 2400A/2400B the modem signal sample rate is 48kHz.  To
listen to the modem tones from FreeDV 2400B, or play them into a FM HT
mic input:
```
$ ./freedv_tx 2400B ../../raw/ve9qrp_10s.raw - | play -t raw -r 48000 -s -2 -
```
Simulate FreeDV 2400B passing through a 300 to 3000 Hz audio path using sox to filter:
```
$  ./freedv_tx 2400B ../../raw/ve9qrp_10s.raw - | sox -t raw -r 48000 -s -2 - -t raw - sinc 300-3000 | ./freedv_rx 2400B - - | play -t raw -r 8000 -s -2 -
```

## FreeDV 2020 support (building with LPCNet)

1. Build codec2 initially without LPCNet
   ```
   $ cd ~
   $ git clone https://github.com/drowe67/codec2.git
   $ cd codec2 && mkdir build_linux && cd build_linux
   $ cmake ../
   $ make
   ```

1. Build LPCNet:
   ```
   $ cd ~
   $ git clone https://github.com/drowe67/LPCNet
   $ cd LPCNet && mkdir build_linux && cd build_linux
   $ cmake -DCODEC2_BUILD_DIR=~/codec2/build_linux ../ 
   $ make
   ```

1. (Re)build Codec 2 with LPCNet support:
   ```
   $ cd ~/codec2/build_linux && rm -Rf *
   $ cmake -DLPCNET_BUILD_DIR=~/LPCNet/build_linux ..
   $ make
   ```

## FreeDV 2020 tests with FreeDV API

```
$ cat ~/LPCNet/wav/wia.wav | ~/LPCNet/build_linux/src/lpcnet_enc -s | ./ofdm_mod --ts 0.0205 --nc 31 --ldpc 2 --verbose 1 -p 312 | ./ofdm_demod --ts 0.0205 --nc 31 --verbose 1 --ldpc 2 -p 312 | ~/LPCNet/build_linux/src/lpcnet_dec -s | aplay -f S16_LE -r 16000
```
Listen the reference tx:
```
$ cat ~/LPCNet/wav/wia.wav | ~/LPCNet/build_linux/src/lpcnet_enc -s | ./ofdm_mod --ts 0.0205 --nc 31 --ldpc 2 --verbose 1 -p 312 | aplay -f S16_LE
```

Listen the freedv_tx:
```
$ ./freedv_tx 2020 ~/LPCNet/wav/wia.wav - | aplay -f S16_LE
```

FreeDV API tx, with reference rx from above:
```
$ ./freedv_tx 2020 ~/LPCNet/wav/wia.wav - | ./ofdm_demod --ts 0.0205 --nc 31 --verbose 1 --ldpc 2 -p 312 | ~/LPCNet/build_linux/src/lpcnet_dec -s | aplay -f S16_LE -r 16000
```

FreeDV API tx and rx:
```
$ ./freedv_tx 2020 ~/LPCNet/wav/all.wav - | ./freedv_rx 2020 - - | aplay -f S16_LE -r 16000
$ ./freedv_tx 2020 ~/LPCNet/wav/all.wav - --testframes | ./freedv_rx 2020 - /dev/null --testframes -vv
```

Simulated HF slow fading channel, 10.8dB SNR:
```
$ ./freedv_tx 2020 ~/LPCNet/wav/all.wav - | ./cohpsk_ch - - -30 --Fs 8000 --slow | ./freedv_rx 2020 - - | aplay -f S16_LE -r 16000
```
It falls down quite a bit with fast fading (--fast):

AWGN (noise but no fading) channel, 2.8dB SNR:
```
$ ./freedv_tx 2020 ~/LPCNet/wav/all.wav - | ./cohpsk_ch - - -22 --Fs 8000 | ./freedv_rx 2020 - - | aplay -f S16_LE -r 16000
```

## TODO

1. Link to this doc from freedv-gui user manual, rowetel/codec2 page, freedv/org
1. README_fdmdv.txt -> .md
1. nice image or two
1. TODO list of relevant files/data demos

