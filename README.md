# Codec 2 README

Codec 2 is an open source (LGPL 2.1) low bit rate speech codec:

http://rowetel.com/codec2.html

Also included:

  + FreeDV API source code.  FreeDV is an open source digital voice
    protocol that integrates the modems, codecs, and FEC
  + FDMDV DPSK modem (README_fdmdv) for HF channels
  + Coherent PSK ((README_cohpsk) for HF channels
  + Non-coherent FSK modem (README_fsk)
  + Coherent OFDM modem for HF channels (README_ofdm)
  + software for High Altitude Balloon image and telemetry reception

## Quickstart

Also see [INSTALL](INSTALL) for more general building and installing instructions. 

1/ Build Codec 2:
```
$ cd codec2
$ mkdir build_linux
$ cd build_linux
$ cmake ..
$ make
```

2/ Listen to Codec 2:
```
$ ./src/c2demo ../raw/hts1a.raw hts1a_c2.raw
$ play -t raw -r 8000 -e signed-integer -b 16 ../raw/hts1a.raw
$ play -t raw -r 8000 -e signed-integer -b 16 ./hts1a_c2.raw
```
3/ Compress, decompress and then play a file:

   using 2400 bps bit rate encoding
```
$ ./src/c2enc 2400 ../raw/hts1a.raw hts1a_c2.bit
$ ./src/c2dec 2400 hts1a_c2.bit hts1a_c2_2400.raw 
```
   which can be played with
```
$ play -t raw -r 8000 -e signed-integer -b 16 ./hts1a_c2_2400.raw
```
   using 700C bps bit rate encoding
```
$ ./src/c2enc 700C ../raw/hts1a.raw hts1a_c2.bit
$ ./src/c2dec 700C hts1a_c2.bit hts1a_c2_700.raw
```
   which can be played with
```
$ play -t raw -r 8000 -e signed-integer -b 16 ./hts1a_c2_700.raw
```
4/ If you prefer a one-liner without saving to files:
```
$ ./src/c2enc 1300 ../raw/hts1a.raw - | ./src/c2dec 1300 - - | play -t raw -r 8000 -b 16 -e signed-integer -
```
   Same at 450 bit/s:
```
$ ./src/c2enc 450 ../raw/ve9qrp.raw - | ./src/c2dec 450 - - | play -t raw -r 8000 -e signed-integer -b 16 -
```
   Please note that 450PWB (pseudo-wideband) can be chosen for decoding, providing a bandwidth extension to 8kHz/16ksps from a 4kHz/8ksps encoded 450bit/s file:
```
$ ./src/c2enc 450 ../raw/ve9qrp.raw - | ./src/c2dec 450PWB - - | play -t raw -r 16000 -e signed-integer -b 16 -
```
## Programs

+ c2demo encodes a file of speech samples, then decodes them and
  saves the result.

+ c2enc encodes a file of speech samples to a compressed file of
  encoded bits.  c2dec decodes a compressed file of bits to a file of
  speech samples.

+ c2sim is a simulation/development version of Codec 2.  It allows
  selective use of the various Codec 2 algorithms.  For example
  switching phase modelling or LSP quantisation on and off.

+ freedv_tx/freedv_rx are command line implementations of the FreeDV
  protocol, which combines Codec 2, modems, and Forward Error
  Correction (FEC).
  
+ cohpsk_* are coherent PSK (COHPSK) HF modem command line programs.

+ drs232, drs232_ldpc, and horus_demod are used for receiving images
  and telemetry from high altitude balloons (Project Horus Wenet,
  Horus Binary protocol)

+ fdmdv_* are differential PSK HF modem command line programs (README_fdmdv).

+ fsk_* are command line programs for a non-coherent FSK modem (README_fsk).

+ ldpc_* are LDPC encoder/decoder command line programs, based on the CML library.

+ ofdm_* are OFDM PSK HF modem command line programs (README_ofdm).

## FreeDV 2020 support (building with LPCNet)

NOTE: Instructions assume you are creating a build_linux directory from within
      the source directory. Adjust paths as needed if this is not the case.

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

## Building and Running Unit Tests

CTest is used as a test frame work, with support from GNU Octave
scripts.

1/ Install GNU Octave on Ubuntu with:
```
$ sudo apt install octave octave-control octave-parallel octave-signal octave-specfun
```
  (see also Octave section below)
  
2/ To build and run the tests:
```
$ cd ~/codec2
$ rm -Rf build_linux && mkdir build_linux
$ cd build_linux
$ cmake -DCMAKE_BUILD_TYPE=Debug ..
$ make all test
```
3/ To just run tests without rebuilding:
```
$ make test
```
4/ To get a verbose run (e.g. for test debugging):
```
$ ctest -V
```
5/ To just run a single test:
```
$ ctest -R test_OFDM_modem_octave_port
```
5/ To list the available tests:
```
$ ctest -N
```
## Directories
```
cmake       - cmake support files
misc        - misc C programs that have been useful in development,
              not reqd for Codec 2 release. Part of Debug build.
octave      - Octave scripts used to support development
script      - shell scripts for playing and converting raw files
src         - C source code for Codec 2, FDMDV modem, COHPSK modem, FreeDV API
raw         - speech files in raw format (16 bits signed linear 8 kHz)
stm32       - STM32F4 microcontroller and SM1000 FreeDV Adaptor support
unittest    - Code to perform and support testing. Part of Debug build.
wav         - speech files in wave file format
```
## GDB and Dump Files

1/ To compile with debug symbols for using gdb:
```
$ cd ~/codec2
$ rm -Rf build_linux && mkdir build_linux
$ cd build_linux
$ CFLAGS=-g cmake ..
$ make
```
2/ For dump file support (dump data from c2sim for input to Octave
development scripts):
```
$ cd ~/codec2
$ rm -Rf build_linux && mkdir build_linux
$ cd build_linux
$ CFLAGS=-DDUMP cmake ..
$ make
```
## Building for Windows on a Linux machine

On Ubuntu 17:
```
$ sudo apt-get install mingw-w64
$ mkdir build_windows && cd build_windows
$ cmake .. -DCMAKE_TOOLCHAIN_FILE=/home/david/freedv-dev/cmake/Toolchain-Ubuntu-mingw32.cmake -DUNITTEST=FALSE -DGENERATE_CODEBOOK=/home/david/codec2/build_linux/src/generate_codebook 
$ make
```
## Building for Windows on a Windows machine

 mkdir build_windows (Or what ever you want to call your build dir)
 cmake -G "MinGW Makefiles" -D CMAKE_MAKE_PROGRAM=mingw32-make.exe
 Or if you use ninja for building cmake -G "Ninja" ..
 mingw32-make or ninja  depends on what you used in the last command
 wait for it to build.

## Octave Packages

To run the Octave scripts the following libraries are required:

Package Name  | Version | Installation directory
--------------|---------|-----------------------
control *     |   2.6.2 | /usr/share/octave/packages/control-2.6.2
general *     |   1.3.4 | /usr/share/octave/packages/general-1.3.4
parallel *    |   2.2.0 | /usr/share/octave/packages/parallel-2.2.0
plot *        |   1.1.0 | /usr/share/octave/packages/plot-1.1.0
signal *      |   1.2.2 | /usr/share/octave/packages/signal-1.2.2
specfun *     |   1.1.0 | /usr/share/octave/packages/specfun-1.1.0

These can be installed using your systems package management system or
the Octave package management system.  The version number of each
package is not important.

On Ubuntu install with:
```
$ sudo apt install octave octave-control octave-parallel octave-signal octave-specfun
```
## FreeDV API

See freedv_api.h and freedv_api.c, and the demo programs freedv_tx &
freedv_rx.  Quickstart demo using FreeDV 1600:
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

