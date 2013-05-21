README.txt
codec2-dev stm32f4 
David Rowe May 2013

Getting Started
-------------------------

Describe what gdb_stdio does, describe what UT does. Where raw files
end up.  Dump files and how to use them.

Install arm toolchain binary

Download and untar STM32F4xx_DSP_StdPeriph_Lib_V1.1.0

build codec2 unit test, describe what it does
make

1. Patching and build stlink:

   TBC

2. start st-util:

   ~/stlink$ sudo ./st-util -f /home/david/codec2-dev/stm32/stm32f4_codec2.elf

3. Start gdb:

   ~/codec2-dev/stm32$ ~/sat/bin/arm-none-eabi-gdb stm32f4_codec2.elf

   (gdb) tar ext :4242

   (gdb) load
   `/home/david/codec2-dev/stm32/fft_test.elf' has changed; re-reading symbols.
    Loading section .isr_vector, size 0x188 lma 0x8000000
    Loading section .text, size 0x1a4b4 lma 0x8000188
    Loading section .data, size 0x28f0 lma 0x801a63c
    Start address 0x800a885, load size 118572
    Transfer rate: 13 KB/sec, 10779 bytes/write.
 
4. Power cycle Discovery.

5. ctrl-C to stop st-util, then restart st-util

6. Back to gdb:

    (gdb) tar ext :4242
    A program is being debugged already.  Kill it? (y or n) y
    Remote connection closed
    (gdb) tar ext :4242
    Remote debugging using :4242
    Reset_Handler () at lib/startup_stm32f4xx.s:69
    69	  movs  r1, #0

   (gdb) c
   Continuing.


Process
-------

1. Profiling macros.

2. enable DUMP to dump files, note proofiling times will be corrupted
by this due to latency in talking to Host

3. Compare outputs using octave/diff_codec.  Worked example:

diff_codec("~/stlink/ref/hts1a_out_1300.raw", "~/stlink/hts1a_out_1300.raw","~/stlink/stm32f4", "~/stlink/ref/stm32f4")


Gotcha
------

using printf rather than gdb_stdio_printf, regular stdio functions are stubbed out so will link, just nothing will happen.

