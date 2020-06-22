# FreeDV README

FreeDV is an open source digital voice protocol that integrates the modems, codecs, and FEC. FreeDV is available as a GUI application, an open source library (FreeDV API), and in hardware (the SM1000 FreeDV adaptor).

## Modes

TODO: Table with name, number of carriers, type of modulation, FEC, typical applications, Codec employed.

## Quickstart

## Programs

+ freedv_tx/freedv_rx are command line implementations of the FreeDV protocol, which combines Codec 2, modems, and Forward Error Correction (FEC).
  
TODO list of relevant files/data demos

## FreeDV 2020 support (building with LPCNet)

NOTE: Instructions assume you are creating a build_linux directory from within the source directory. Adjust paths as needed if this is not the case.

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

### FreeDV 2020 tests with FreeDV API

Reference: Plugging together lpcnet_enc -> ofdm_mod -> ofdm_demod -> lpcnet_dec:
```
$ cat ~/LPCNet/wav/wia.wav | ~/LPCNet/build_linux/src/lpcnet_enc -s | ./ofdm_mod --ts 0.0205 --nc 31 --ldpc 2 --verbose 1 -p 312 | ./ofdm_demod --ts 0.0205 --nc 31 --verbose 1 --ldpc 2 -p 312 | ~/LPCNet/build_linux/src/lpcnet_dec -s | aplay -f S16_LE -r 16000
```
We are trying to integrate all of the above into FreeDV API.

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
It falls down quite a bit with fast fading (--fast).  We'll work on that.

AWGN (noise but no fading) channel, 2.8dB SNR:
```
$ ./freedv_tx 2020 ~/LPCNet/wav/all.wav - | ./cohpsk_ch - - -22 --Fs 8000 | ./freedv_rx 2020 - - | aplay -f S16_LE -r 16000
```

## FreeDV API

See freedv_api.h and freedv_api.c, and the demo programs freedv_tx & freedv_rx.  Quickstart demo using FreeDV 1600:
```
$ ./freedv_tx 1600 ../../raw/hts1.raw - | ./freedv_rx 1600 - - | play -t raw -r 8000 -s -2 -q -
$ cat freedv_rx_log.txt
```
## FreeDV 2400A and 2400B modes

FreeDV 2400A and FreeDV 2400B are modes designed for VHF radio.
FreeDV 2400A is designed for SDR radios (it has a 5 kHz RF bandwidth),
however FreeDV 2400B is designed to pass through commodity FM radios.

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

## Links:

+ FreeDV 2400A blog post ...: http://www.rowetel.com/?p=5119
+ FreeDV 2400A & 2400B demos: http://www.rowetel.com/?p=5219

