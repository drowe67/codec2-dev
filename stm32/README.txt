README.txt
codec2-dev stm32f4 
David Rowe May 2013

TODO
 + Describe what gdb_stdio does, describe what UT does.
 + Where raw files end up.  
 + Dump files and how to use them.
 + check if "CFLAGS:  -mlittle-endian -mthumb -mthumb-interwork" needed

Getting Started
-------------------------

. Install arm toolchain binary

   $ cd ~
   $ wget https://launchpadlibrarian.net/126639661/gcc-arm-none-eabi-4_7-2012q4-20121208-linux.tar.bz2
   $ tar xjf gcc-arm-none-eabi-4_7-2012q4-20121208-linux.tar.bz2

. Build codec2 unit test:

   $ cd codec2_dec/stm
   In the Makefile edit the BINPATH variable for your toolchain location
   $ make

. Patching and build stlink:

   $ cd ~
   $ git clone https://github.com/texane/stlink.git
   $ cd stlink
   ~/stlink$ git checkout bbecbc1e81b15b85829149424d048d96bd844939
   ~/stlink$ patch -p0 < ~/codec2-dev/stm32/stlink/stlink.patch
   ~/stlink$ cp ~/codec2-dev/stm32/stlink/elfsym.* .
   ~/stlink$ ./autogen.sh
   ~/stlink$ ./configure
   ~/stlink$ make 

. Place a copy of hts1a.raw in the stlink directory and start st-util:

   ~/stlink$ cp ~/codec2-dev/raw/hts1a .
   ~/stlink$ sudo ./st-util -f /home/david/codec2-dev/stm32/stm32f4_codec2.elf

. In _another_ console start gdb:

   $ ~/codec2-dev/stm32$ ~/sat/bin/arm-none-eabi-gdb stm32f4_codec2.elf

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
  
  ~/stlink$ play -r 8000 -s -2 hts1a_out.raw

Process
-------

1. Profiling macros.

2. enable DUMP variable in Makefile to dump files, note profiling
times will be corrupted by this due to latency in talking to Host

3. Compare outputs using octave/diff_codec.  Worked example:

diff_codec("~/stlink/ref/hts1a_out_1300.raw", "~/stlink/hts1a_out_1300.raw","~/stlink/stm32f4", "~/stlink/ref/stm32f4")


Gotcha
------

Using printf rather than gdb_stdio_printf, regular stdio functions are stubbed out so will link, just nothing will happen.

