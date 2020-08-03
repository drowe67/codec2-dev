# FreeDV Technology

FreeDV is an open source digital voice protocol that integrates modems, codecs, and FEC.

On transmit, FreeDV converts speech to a modem signal you can send over a radio channel.  On receive, FreeDV takes off air modem signals and converts them to speech samples.

FreeDV is available as a GUI application, an open source library (FreeDV API), and in hardware (the SM1000 FreeDV adaptor).  FreeDV is part of the Codec 2 project.

This document gives an overview of the technology inside FreeDV, and some additional notes on building/using the FreeDV 2020 and 2400A/2400B modes.

![FreeDV mode knob](http://www.rowetel.com/images/codec2/mode_dv.jpg)

## Reading Further

1. [FreeDV web site](http://freedv.org)
1. [FreeDV GUI User Manual](https://github.com/drowe67/freedv-gui/blob/master/USER_MANUAL.md)
1. [Codec 2](http://rowetel.com/codec2.html)
1. FreeDV can also be used for data [README_data](https://github.com/drowe67/codec2/blob/master/README_data.md)
1. [FreeDV 1600 specification](https://freedv.org/freedv-specification)
1. [FreeDV 700C blog post](http://www.rowetel.com/wordpress/?p=5456)
1. [FreeDV 700D Released blog post](http://www.rowetel.com/wordpress/?p=6103)
1. [FreeDV 2020 blog post](http://www.rowetel.com/wordpress/?p=6747)
1. [FreeDV 2400A blog post](http://www.rowetel.com/?p=5119)
1. [FreeDV 2400A & 2400B](http://www.rowetel.com/?p=5219)
1. Technical information on various modem waveforms in the [modem codec frame design spreadsheet](https://github.com/drowe67/codec2/blob/master/doc/modem_codec_frame_design.ods)
1. [Modems for HF Digital Voice Part 1](http://www.rowetel.com/wordpress/?p=5420)
1. [Modems for HF Digital Voice Part 2](http://www.rowetel.com/wordpress/?p=5448)
1. [FDMDV modem README](README_fdmdv.md)
1. [OFDM modem README](README_ofdm.md)
1. Many blog posts in the [rowetel.com blog archives](http://www.rowetel.com/?page_id=6172)

## FreeDV HF Modes

These are designed for use with a HF SSB radio.

| Mode | Date | Codec | Modem | RF BW | Raw bits/s | FEC | Text | SNR min | Multipath |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1600 | 2012 | Codec2 1300 | DQPSK + pilot | 1125 | 1600 | Golay (23,12) | Y | 4 | poor |
| 700C | 2017 | Codec2 700C | coherent QPSK + diversity | 1500 | 1400 | - | - | 2 | good |
| 700D | 2018 | Codec2 700C | coherent OFDM/QPSK | 1000 | 1900 | LDPC (224,112)  | Y | -2 | fair |
| 2020 | 2019 | LPCNet 1733 | coherent OFDM/QPSK | 1600 | 3000 | LDPC (504,396)  | Y | 2 | poor |

Notes:

1. *Raw bits/s* is the number of payload bits/s carried over the channel by the modem.  This consists of codec frames, FEC parity bits, unprotected text, and synchronisation information such as pilot and unique word bits.  The estimates are open to interpretation for the OFDM waveforms due to pilot symbol and cyclic prefix considerations (see spreadsheet).

1. *RF BW* is the bandwidth of the RF signal over the air.  FreeDV is more bandwidth efficient than SSB.

1. *Multipath* is the relative resilience of the mode to multipath fading, the biggest problem digital voice faces on HF radio channels.  Analog SSB would be rated as "good".

1. *Text* is a side channel for low bit rate text such as your location and call sign.  It is generally unprotected by FEC, and encoded with varicode.

1. *SNR Min* is for an AWGN channel (no multipath/fading).

1. All of the modems use multiple parallel carriers running at a low symbol rate of around 50 Hz.  This helps combat the effects of multipath channels.

1. Some of the Codec 2 modes (2400/1300/700C etc) happen to match the name of a FreeDV mode.  For example FreeDV 700C uses Codec 2 700C for voice compression. However FreeDV 700D *also* uses Codec 2 700C for voice compression, but has a very different modem waveform to FreeDV 700C.  Sorry for the confusing nomenclature.

1. Coherent demodulation gives much better performance than differential, at the cost of some additional complexity.  Pilot symbols are transmitted regularly to allow the demod to estimate the reference phase of each carrier.

1. The 1600 and 700C waveforms use parallel tone modems, later modes use OFDM.  OFDM gives tighter carrier packing which allows higher bit rates, but tends to suffer more from frequency offsets and delay spread.

1. At medium to high SNRs, FreeDV 700C performs well (better than 700D) on fast fading multipath channels with large delay spread due its parallel tone design and high pilot symbol rate.  It employs transmit diversity which delivers BER performance similar to modes using FEC.  FreeDV 700C also has a short frame (40ms), so syncs fast with low latency.  Fast sync is useful on marginal channels that move between unusable and barely usable.

1. FreeDV 700D uses an OFDM modem and was optimised for low SNR channels, with strong FEC but a low pilot symbol rate and modest (2ms) cyclic prefix which means its performance degrades on multipath channels with fast (> 1Hz) fading.  The use of strong FEC makes this mode quite robust to other channel impairments, such as static crashes, urban HF noise, and in-band interference.

1. FEC was added fairly recently to FreeDV modes.  The voice codecs we use work OK at bit error rates of a few %, and packet error rates of 10%. Raw bit error rates on multipath channels often exceed 10%.  For reasonable latency (say 40ms) we need small codewords. Thus to be useful we require a FEC code that works at over 10% raw BER, has 1% output (coded) bit error rate, and a codeword of around 100 bits.  Digital voice has unusual requirements, most FEC codes are designed for data which is intolerant of any bit errors, and few operate over 10% raw BER.  Powerful FEC codes have long block lengths (1000's of bits) which leads to long latency.  However LDPC codes come close, and can also "clean up" other channel errors caused by static and interference.  The use of OFDM means we now have "room" for the extra bits required for FEC, so there is little cost in adding it, apart from latency.

## FreeDV VHF Modes

These modes use constant amplitude modulation like FSK or FM, and are designed for VHF and above.  However 800XA can be run over HF or VHF on a SSB radio.

| Mode | Date | Codec2 | Modem | RF BW | Raw bits/s | FEC | Text |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 2400A | 2016 | 1300 | 4FSK | 2400 | 5kHz | Golay (23,12) | Y |
| 2400B | 2016 | 1300 | baseband/analog FM | analog FM | 2400 | Golay (23,12) | Y |
| 800XA | 2017 | 700C |  4FSK | 2000 | 800 | - | N |

## FreeDV API

See [freedv_api.h](src/freedv_api.h) and [freedv_api.c](src/freedv_api.c), and the command line demo programs [freedv_tx.c](src/freedv_tx.c) & [freedv_rx.c](src/freedv_rx.c).  Quickstart demo using FreeDV 1600:
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
$ cat ~/LPCNet/wav/wia.wav | ~/LPCNet/build_linux/src/lpcnet_enc -s | ./ofdm_mod --mode 2020 --ldpc --verbose 1 -p 312 | ./ofdm_demod --mode 2020 --verbose 1 --ldpc -p 312 | ~/LPCNet/build_linux/src/lpcnet_dec -s | aplay -f S16_LE -r 16000
```
Listen the reference tx:
```
$ cat ~/LPCNet/wav/wia.wav | ~/LPCNet/build_linux/src/lpcnet_enc -s | ./ofdm_mod --mode 2020 --ldpc --verbose 1 -p 312 | aplay -f S16_LE
```

Listen the freedv_tx:
```
$ ./freedv_tx 2020 ~/LPCNet/wav/wia.wav - | aplay -f S16_LE
```

FreeDV API tx, with reference rx from above:
```
$ ./freedv_tx 2020 ~/LPCNet/wav/wia.wav - | ./ofdm_demod --mode 2020 --verbose 1 --ldpc -p 312 | ~/LPCNet/build_linux/src/lpcnet_dec -s | aplay -f S16_LE -r 16000
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

## Further work

1. ~~Spell check/proofread~~
1. Screen shots of each modes (waterfall and spectrum), discussing aspects of waveforms
1. ~~Link to this doc from freedv-gui user manual, rowetel/codec2 page, freedv.org~~
1. ~~README_fdmdv.txt -> .md~~
1. nice image or two
1. table of source files
