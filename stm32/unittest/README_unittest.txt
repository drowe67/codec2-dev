README_unittest.txt
Don Reid 2018/2019

This is the unittest system for the stm32 implementation of codec2/FreeDV
It is currently only working on linux systems. It requires a STM32F4xx processor 
development board connected to/having a ST-LINK , e.g. 
a STM32F4 Discovery board. 

Quickstart
==========

   $ cd <path_to_codec2>/stm32/unittest
   $ ./scripts/run_all_stm32_tests


Objectives
==========

   It is important to have a robust set of tests for the functionality
   of any software project, especially when there are multiple developers.

   It is easier to test small units of code first as stand alone functions.

   The more automated the test system is the easier it is to run and thus
   the more likely people are to run it.


Directory Structure
===================

   scripts     Where scripts for this unittest system are found and run

   src         Where the stm32 sources and binaries are, Makefile to build tests

   lib         Where includable files for this unittest system are found
      /python     python library files
      /octave     octave library files

   test_run    Each test is run in a subdirectory here.


Debug and semihosting
=====================

   These tests use a recent version of the st-util program.

   The source can be downloaded from:

       https://github.com/texane/stlink

   After compiling it can be installed anywhere.  The program is in
   build/Release/src/gdbserver/st-util.

   The newlib stdio functions (open, fread, fwrite, flush, fclose, etc.) send
   some requests that this tool does not recognize and those messages will appear
   in the output of st-util.  They can be ignored.


Building and Running the stm32 Unit Tests
=========================================

0/ Build codec2 for linux and stm32 

   Run the linux cmake build using /path/to/codec2/build_linux 
   as build directory, see instructions in the codec2 top level directory.

   Run the stm32 cmake build using /path/to/codec2/stm32/build_stm32 
   as build directory, see instructions in the stm32 directory.
   
1/ Build stlink:

  $ cd ~
  $ git clone https://github.com/texane/stlink
  $ cd stlink
  $ make
  $ sudo cp ./etc/udev/rules.d/49-stlinkv2.rules /etc/udev/rules.d/
  $ sudo udevadm control --reload-rules

1a/ Add the st-util util to your $PATH

1b/ Plug in a stm32 development board and test:

  $ st-util

  st-util 1.4.0-47-gae717b9
  2018-12-29T06:52:16 INFO usb.c: -- exit_dfu_mode
  2018-12-29T06:52:16 INFO common.c: Loading device parameters....
  2018-12-29T06:52:16 INFO common.c: Device connected is: F4 device, id 0x10016413
  2018-12-29T06:52:16 INFO common.c: SRAM size: 0x30000 bytes (192 KiB), Flash: 0x100000 bytes (1024 KiB) in pages of 16384 bytes
  2018-12-29T06:52:16 INFO gdb-server.c: Chip ID is 00000413, Core ID is  2ba01477.
  2018-12-29T06:52:16 INFO gdb-server.c: Listening at *:4242...

2/ To run a single test:

   $ cd codec2-dev/stm32/unittest
   $ ./scripts/run_stm32_test <name_of_test> <test_option>

for example 

   $ ./scripts/run_stm32_tst tst_ofdm_demod quick

3/ To run a test set (codec2, ofdm, ldpc):

   $ cd codec2-dev/stm32/unittest
   $ ./scripts/run_all_<set_name>_tests

for example 

   $ ./scripts/run_all_ldpc_tests

4/ To run ALL tests, see "Quickstart" above


# vi:set ts=3 et sts=3:
