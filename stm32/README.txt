README.txt
codec2 support for the stm32f4
David Rowe Jan 2016

Quickstart
==========

The Makefile generates several programs used in Codec 2 development on
the STM32F4, including sm1000.bin, the firmware for the SM1000.

1. Install the toolchain, on Ubuntu 14 this is:

   $ sudo apt-get install gcc-arm-none-eabi

   On Fedora:

   $ sudo dnf install stlink gcc-arm-linux-gnu


2. $ make (and cross your fingers)

Not quite so Quickstart
=======================

Note: This section needs some editing.  It deals with the running the
semi-hosting used for development system.

gdb_stdio system
----------------

stutil contains a gdb server that talks to the target firmware.
stutil has been patched to allow "semihosting": stdio requests on the
target are re-directed to the host PC.  So if you call printf on the
target, it appears on the host PC console.  With printf/fread/fwrite
and gdb it makes developing on bare metal just like developing on any
command line gcc system.

The root path for files accessed by the target is the path st-util is
run from.

Getting Started
---------------

. Install arm toolchain binary

   $ cd ~
   $ wget https://launchpad.net/gcc-arm-embedded/4.7/4.7-2013-q1-update/+download/gcc-arm-none-eabi-4_7-2013q1-20130313-linux.tar.bz2
   $ tar xjf gcc-arm-none-eabi-4_7-2013q1-20130313-linux.tar.bz2

. Build codec2 unit test:

   $ cd codec2_dev/stm32
   In Makefile: edit the BINPATH variable for your toolchain location
                edit PERIPHLIBVER for the current version of the peripheral
                 library, currently V1.3.0
                delete power_ut.elf target from the all: target
   $ make
   after make downloads the peripheral library, you will get a compile error:
   #error "Please select first the target STM32F4xx device used in your application (in stm32f4xx.h file)"
   edit STM32F4xx_DSP_StdPeriph_Lib_V1.3.0/Libraries/CMSIS/Device/ST/STM32F4xx/Include/stm32f4xx.h
   and at the top of that file, uncomment the appropriate line for your target
   processor.
   $ make

. Build stlink:

   $ cd ~
   $ git clone https://github.com/shenki/stlink.git
   $ cd stlink
   ~/stlink$ sudo apt-get install libusb-1.0-0-dev libelf-dev automake
   ~/stlink$ ./autogen.sh
   ~/stlink$ ./configure
   ~/stlink$ make

. Place a copy of hts1a.raw in the stlink directory and start st-util:

   ~/stlink$ cp ~/codec2-dev/raw/hts1a.raw stm_in.raw
   ~/stlink$ sudo ./st-util -f /home/david/codec2-dev/stm32/stm32f4_codec2.elf

. In _another_ console start gdb:

   $ ~/codec2-dev/stm32$ ~/gcc-arm-none-eabi-4_7-2013q1/bin/arm-none-eabi-gdb stm32f4_codec2.elf

   (gdb) tar ext :4242

   (gdb) load
   `/home/david/codec2-dev/stm32/fft_test.elf' has changed; re-reading symbols.
    Loading section .isr_vector, size 0x188 lma 0x8000000
    Loading section .text, size 0x1a4b4 lma 0x8000188
    Loading section .data, size 0x28f0 lma 0x801a63c
    Start address 0x800a885, load size 118572
    Transfer rate: 13 KB/sec, 10779 bytes/write.

. Power cycle Discovery.

. Stop st-util using ctrl-C, then restart st-util

. Back to gdb:

    (gdb) tar ext :4242
    A program is being debugged already.  Kill it? (y or n) y
    Remote connection closed
    (gdb) tar ext :4242
    Remote debugging using :4242
    Reset_Handler () at lib/startup_stm32f4xx.s:69
    69	  movs  r1, #0

   (gdb) c
   Continuing.

. gdb will prints various debug messages, and the codec output file
  will be written to ~/stlink.

  ~/stlink$ play -r 8000 -s -2 stm_out.raw

Process
-------

1. Profiling macros, grep on TIMER_xxxx

2. Enable DUMP variable in Makefile to dump files, note profiling
   times will be corrupted by this due to latency in talking to Host

3. Compare outputs using octave/diff_codec.  Example:

   octave:> diff_codec("~/stlink/ref/hts1a_out_1300.raw", "~/stlink/hts1a_out_1300.raw","~/stlink/stm32f4", "~/stlink/ref/stm32f4")

Gotcha
------

Using printf rather than gdb_stdio_printf, regular stdio functions are
stubbed out so will link, just nothing will happen.

TODO
----

 + check if "CFLAGS:  -mlittle-endian -mthumb -mthumb-interwork" needed
 + double check if _fini hack is OK (src/init.c)
