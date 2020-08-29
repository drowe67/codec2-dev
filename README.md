# Codec 2 README

Codec 2 is an open source (LGPL 2.1) low bit rate speech codec: http://rowetel.com/codec2.html

Also included:

  + The FreeDV API for digital voice over radio. FreeDV is an open source digital voice protocol that integrates modems, codecs, and FEC [README_freedv](README_freedv.md)
  + APIs for raw and Ethernet packet data over radio [README_data](README_data.md)
  + High performance coherent OFDM modem for HF channels [README_ofdm](README_ofdm.md)
  + High performance non-coherent FSK modem [README_fsk](README_fsk.md)
  + An STM32 embedded version of FreeDV 1600/700D for the [SM1000](stm32/README.md)
  + Coherent PSK modem [README_cohpsk](README_cohpsk.md) for HF channels
  + FDMDV DPSK modem [README_fdmdv](README_fdmdv.md) for HF channels

## Quickstart

1. Install packages (Debian/Ubuntu):
   ```
   $ sudo apt install git build-essential cmake
   ```
   Fedora/RH distros:
   ```
   $ sudo dnf groupinstall "Development Tools" "C Development Tools and Libraries"
   $ sudo dnf install cmake
   ```
   
1. Build Codec 2:
   ```
   $ git clone git@github.com:drowe67/codec2.git
   $ cd codec2
   $ mkdir build_linux
   $ cd build_linux
   $ cmake ..
   $ make
   ```

1. Listen to Codec 2 (run from ```codec2/build_linux```):
   ```
   $ ./src/c2demo ../raw/hts1a.raw hts1a_c2.raw
   $ aplay -f S16_LE ../raw/hts1a.raw
   $ aplay -f S16_LE hts1a_c2.raw
   ```
1. Compress, decompress and then play a file using Codec 2 at 2400 bit/s:
   ```
   $ ./src/c2enc 2400 ../raw/hts1a.raw hts1a_c2.bit
   $ ./src/c2dec 2400 hts1a_c2.bit hts1a_c2_2400.raw 
   ```
   which can be played with:
   ```
   $ aplay -f S16_LE hts1a_c2_2400.raw
   ```
   Or using Codec 2 using 700C (700 bits/s):
   ```
   $ ./src/c2enc 700C ../raw/hts1a.raw hts1a_c2.bit
   $ ./src/c2dec 700C hts1a_c2.bit hts1a_c2_700.raw
   $ aplay -f S16_LE hts1a_c2_700.raw
   ```
1. If you prefer a one-liner without saving to files:
   ```
   $ ./src/c2enc 1300 ../raw/hts1a.raw - | ./src/c2dec 1300 - - | aplay -f S16_LE
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

## Programs

+ ```c2demo``` encodes a file of speech samples, then decodes them and saves the result.

+ ```c2enc``` encodes a file of speech samples to a compressed file of encoded bits.  ```c2dec``` decodes a compressed file of bits to a file of speech samples.

+ ```c2sim``` is a simulation/development version of Codec 2.  It allows selective use of the various Codec 2 algorithms.  For example switching phase modelling or quantisation on and off.

+ ```freedv_tx``` & ```freedv_rx``` are command line implementations of the FreeDV protocol, which combines Codec 2, modems, and Forward Error Correction (FEC).
  
+ ```cohpsk_*``` are coherent PSK (COHPSK) HF modem command line programs.

+ ```fdmdv_*``` are differential PSK HF modem command line programs (README_fdmdv).

+ ```fsk_*``` are command line programs for a non-coherent FSK modem (README_fsk).

+ ```ldpc_*``` are LDPC encoder/decoder command line programs, based on the CML library.

+ ```ofdm_*``` are OFDM PSK HF modem command line programs (README_ofdm).

## Building and Running Unit Tests

CTest is used as a test framework, with support from [GNU Octave](https://www.gnu.org/software/octave/) scripts.

1. Install GNU Octave and libraries on Ubuntu with:
   ```
   $ sudo apt install octave octave-common octave-signal liboctave-dev gnuplot python3-numpy sox valgrind
   ```
1. Install CML library with instructions at the top of [```octave/ldpc.m```](octave/ldpc.m)

1. To build and run the tests:
   ```
   $ cd ~/codec2
   $ rm -Rf build_linux && mkdir build_linux
   $ cd build_linux
   $ cmake -DCMAKE_BUILD_TYPE=Debug ..
   $ make all test
   ```

1. To just run tests without rebuilding:
   ```
   $ ctest
   ```

1. To get a verbose run (e.g. for test debugging):
   ```
   $ ctest -V
   ```

1. To just run a single test:
   ```
   $ ctest -R test_OFDM_modem_octave_port
   ```

1. To list the available tests:
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

1. To compile with debug symbols for using gdb:
   ```
   $ cd ~/codec2
   $ rm -Rf build_linux && mkdir build_linux
   $ cd build_linux
   $ CFLAGS=-g cmake ..
   $ make
   ```

1. For dump file support (dump data from c2sim for input to Octave development scripts):
   ```
   $ cd ~/codec2
   $ rm -Rf build_linux && mkdir build_linux
   $ cd build_linux
   $ CFLAGS=-DDUMP cmake ..
   $ make
   ```

## Building for Windows on a Linux machine

On Ubuntu 17 and above:
   ```
   $ sudo apt-get install mingw-w64
   $ mkdir build_windows && cd build_windows
   $ cmake .. -DCMAKE_TOOLCHAIN_FILE=/home/david/freedv-dev/cmake/Toolchain-Ubuntu-mingw32.cmake -DUNITTEST=FALSE -DGENERATE_CODEBOOK=/home/david/codec2/build_linux/src/generate_codebook 
   $ make
   ```
   
## Building for Windows on a Windows machine

To build with MinGW gcc:
 ```
 mkdir build_windows (Or what ever you want to call your build dir)
 cd build_windows
 cmake .. -G "MinGW Makefiles" -D CMAKE_MAKE_PROGRAM=mingw32-make.exe
 Or if you use ninja for building cmake -G "Ninja" ..
 mingw32-make or ninja  depends on what you used in the last command
 wait for it to build.
 ```

To build with Microsoft Visual Studio 2017 or 2019:
 ```
 mkdir build_windows (Or what ever you want to call your build dir)
 cd build_windows
 cmake .. -G "Visual Studio 15 2017 Win64" -DCMAKE_BUILD_TYPE=Release
 or
 cmake .. -G "Visual Studio 16 2019" -A x64 -DCMAKE_BUILD_TYPE=Release
 cmake --build . --config Release --target codec2
 ```
Currently this will just build ./src/Release/codec2.dll, but not the demo apps.
