Codec 2 README
--------------

Codec 2 is an open source 2400 bit/s speech codec (LGPL licensed).
For more information please see:

    http://rowetel.com/codec2.html

Quickstart
----------

1/ Listen to Codec 2:

   $ cd codec2/src
   $ make
   $ ./c2demo ../raw/hts1a.raw hts1a_c2.raw
   $ ../script/menu.sh ../raw/hts1a.raw hts1a_c2.raw

   NOTE: For playback testing, menu.sh requires either the 'play',
   'aplay' or 'ossplay' programs to be installed (see
   http://sox.sourceforge.net/, http://www.alsa-project.org/, or
   http://www.opensound.com/ respectively).

2/ Compress and Decompress a file:

   $ ./c2enc ../raw/hts1a.raw hts1a_c2.bit
   $ ./c2dec hts1a_c2.bit hts1a_c2.raw 

Programs
--------
 
1/ c2demo encodes a file of speech samples, then decodes them and
saves the result.

2/ c2enc encodes a file of speech samples to a compressed file of
encoded bits.

3/ c2dec decodes a compressed file of bits to a file of speech
samples.

4/ c2sim is a simulation/development version of Codec 2.  It allows
selective use of the various Codec 2 algorithms.  For example
switching phase modelling or LSP quantisation on and off.

Debugging
---------

1/ For dump file support:

  $ cd codec2
  $ CFLAGS=-DDUMP ./configure
  $ make clean && make

2/ To use gdb:

  $ $ libtool --mode=execute gdb c2sim

Directories
-----------

  script   - shell scripts for playing and converting raw files
  src      - C source code
  octave   - Octave scripts used for visualising internal signals 
             during development
  raw      - speech files in raw format (16 bits signed linear 8 kHz)
  unittest - unit test source code
  voicing  - hand-estimated voicing files, used for development
  wav      - speech files in wave file format

