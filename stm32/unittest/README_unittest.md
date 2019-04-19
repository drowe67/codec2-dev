# README Codec2 STM32 Unit Test
Don Reid 2018/2019

This is the unittest system for the stm32 implementation of codec2/FreeDV
It is currently only working on linux systems. It requires a STM32F4xx processor 
development board connected to/having a ST-LINK , e.g. 
a STM32F4 Discovery board. 

## Quickstart

* You must have build both build_linux and stm32/build according to the instructions
* You must have numpy for Python3 installed
* You must have either st-util or openocd installed and in you PATH. 
* You must have an arm-none-eabi-gdb install and in the path.

Then do:

```Bash
   cd <path_to_codec2>/stm32/build_stm32
   make test 
   #or
   ctest
```

`ctest` / `make test` will run all tests. `ctest -V`  will provide more output.

You should see a couple of tests executing (and passing). 

Please note that the tests tst_ofdm_demod_AWGN, tst_ofdm_demod_fade, 
tst_ofdm_demod_ldpc_AWGN, tst_ofdm_demod_ldpc_fade fail on certain OS such 
as Ubuntu 18.04.  All other tests MUST pass.


## Objectives

It is important to have a robust set of tests for the functionality
of any software project, especially when there are multiple developers.

It is easier to test small units of code first as stand alone functions.

The more automated the test system is the easier it is to run and thus
the more likely people are to run it.


## Directory Structure

   |Path        |        Description|
   |------------|-------------------|
   |`scripts`   |   Where scripts for this unittest system are found and run
   |`src`       | Where the stm32 sources are
   |`lib`       | Where includable files for this unittest system are found
   |`lib/python`|    python library files
   |`lib/octave`|    octave library files
   |`test_run`  |Each test is run in a subdirectory here.


## Debug and Semihosting

These tests required a connection from the arm-none-eabi-gdb debugger to the
stm32 hardware. For this we use a recent version of OpenOCD or alternatively 
the st-util program.

### OpenOCD
We recommend to use openocd instead of stlink. 

Most linux distributions use too old packages for openocd, so you in general
should build it from the git source. If your test runs fail with error 
messages rearding SYS_FLEN not support (see openocd_stderr.log in the test_run 
directories), your openocd is too old. Make sure to have the right openocd
first in the PATH if you have multiple openocd installed!

It is recommended to build OpenOCD from sources, see below.

### st-util
Most distributions don't have stutil included. Easiest way is to build it from
the github sources.

The source can be downloaded from:

       (https://github.com/texane/stlink)

After compiling it can be installed anywhere, as long as it is in the PATH.  The program is in
`build/Release/src/gdbserver/st-util`.

The newlib stdio functions (open, fread, fwrite, flush, fclose, etc.) send
some requests that this tool does not recognize and those messages will appear
in the output of st-util.  They can be ignored.



## Preparing the Unittest environment

### Build codec2 for linux and stm32 
It is important to build in the plaths defined below, since the test system
assumes some binaries to be located there.

Run the linux cmake build using /path/to/codec2/build_linux 
as build directory, see instructions in the codec2 top level directory.

Run the stm32 cmake build using /path/to/codec2/stm32/build_stm32 
as build directory, see instructions in the stm32 directory.
   
### Prepare stlink (alternative to the recommended OpenOCD, see below):

1. Build from github

```Bash
  git clone https://github.com/texane/stlink
  cd stlink
  make
  sudo cp ./etc/udev/rules.d/49-stlinkv2.rules /etc/udev/rules.d/
  sudo udevadm control --reload-rules
```
3. Add the st-util util to your $PATH, if not installed in the default location ( /usr/local/bin )

4. Plug in a stm32 development board and test:

```
  $ st-util

  st-util 1.4.0-47-gae717b9
  2018-12-29T06:52:16 INFO usb.c: -- exit_dfu_mode
  2018-12-29T06:52:16 INFO common.c: Loading device parameters....
  2018-12-29T06:52:16 INFO common.c: Device connected is: F4 device, id 0x10016413
  2018-12-29T06:52:16 INFO common.c: SRAM size: 0x30000 bytes (192 KiB), Flash: 0x100000 bytes (1024 KiB) in pages of 16384 bytes
  2018-12-29T06:52:16 INFO gdb-server.c: Chip ID is 00000413, Core ID is  2ba01477.
  2018-12-29T06:52:16 INFO gdb-server.c: Listening at *:4242...
```

### Build OpenOCD (recommended over stlink):
OpenOCD needs to be built from the source. If you have successfully build
the linux codec2 binaries, everything required to build OpenOCD is already installed.

1.
    The executable is placed in /usr/local/bin ! Make sure to have no
    other openocd installed (check output of `which openocd` to be 
    /usr/local/bin)

	```Bash
      git clone https://git.code.sf.net/p/openocd/code openocd-code
      cd openocd-code
      ./bootstrap
      ./configure
      sudo make install
      which openocd
	```

2. Plug in a stm32 development board and test:

```
   $ openocd -f board/stm32f4discovery.cfg

   Open On-Chip Debugger 0.10.0+dev-00796-ga4ac5615 (2019-04-12-21:58)
   Licensed under GNU GPL v2
   For bug reports, read
   http://openocd.org/doc/doxygen/bugs.html
   Info : The selected transport took over low-level target control. The results might differ compared to plain JTAG/SWD
   adapter speed: 2000 kHz
   adapter_nsrst_delay: 100
   none separate
   srst_only separate srst_nogate srst_open_drain connect_deassert_srst
   Info : Listening on port 6666 for tcl connections
   Info : Listening on port 4444 for telnet connections
   Info : clock speed 2000 kHz
   Info : STLINK V2J33S0 (API v2) VID:PID 0483:3748
   Info : Target voltage: 2.871855
   Info : stm32f4x.cpu: hardware has 6 breakpoints, 4 watchpoints
   Info : Listening on port 3333 for gdb connections
```

   To run the tests with openocd instead of the default st-util, add `--openocd` to the command lines shown below.

### Install numpy for Python3
Some test are in fact python3 scripts and require the numpy package to be installed,
otherwise some tests will fail.

On Ubuntu:
   ```Bash
   sudo apt-get install python3-numpy 
   ```
 

## Running the stm32 Unit Tests
1. Build Linux and STM32 binaries in codec2/build_linux and codec2/stm32/build_stm32 respectively
   If using OpenOCD pass -DUT_PARAMS=--openocd to your STM32 cmake call!
   
2. Tests can be run using the ctest utility (part  of cmake)
   ```
   $ cd /path/to/codec2/stm32/build_stm32
   $ ctest 
   ```  
   You can pass -V to see more output:
   ```
   $ cd /path/to/codec2/stm32/build_stm32
   $ ctest -V
   ```` 

   You can pass -R <pattern> to run test matching <pattern>. Please note,
   that some test have dependencies and will have to run other tests before
   being executed
   ```
   $ cd /path/to/codec2/stm32/build_stm32
   $ ctest -R ofdm
   ```
   
3. To simply run all tests the `make test` target exists (part  of cmake)
   ```
   $ cd /path/to/codec2/stm32/build_stm32
   $ make test 
   ```
   
4. To run a single test directly:
   ```
   $ cd codec2/stm32/unittest
   $ ./scripts/run_stm32_test <name_of_test> <test_option> --load
   ```
   for example 
   ```
   $ ./scripts/run_stm32_tst tst_ofdm_demod quick --load
   ```
   (Note when running a single test you can choose not to reload the flash every
   time if using the same bits.  The *_all_* scripts manage this themselves.)

5. To run a test set (codec2, ofdm, ldpc):
   ```
   $ cd codec2/stm32/unittest
   $ ./scripts/run_all_<set_name>_tests
   ```
   for example 
   ```
   $ ./scripts/run_all_ldpc_tests
   ```
   
6. To run ALL tests, see "Quickstart" above


# vi:set ts=3 et sts=3:
