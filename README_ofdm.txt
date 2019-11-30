
README_ofdm.txt
David Rowe
Created Mar 2018

Introduction
------------

An Orthogonal Frequency Division Multiplexed (OFDM) modem designed for
digital voice over HF SSB.  Typical configuration for FreeDV 700D is
700 bit/s voice, a rate 0.5 LDPC code, and 1400 bit/s raw data rate
over the channel.

The OFDM modem was first implemented in GNU Octave, then ported to C.
Algorithm development is generally easier in Octave, but for real time
work we need the C version.  Automated units tests ensure the
operation of the Octave and C versions are identical.

Refs
-----

[1] Spreadsheet (OFDM tab) describing the waveform design,
    https://svn.code.sf.net/p/freetel/code/codec2-dev/octave/modem_codec_frame_design.ods

[2] Towards FreeDV 700D, https://www.rowetel.com/?p=5573

[3] FreeDV 700D - First Over The Air Tests, https://www.rowetel.com/?p=5630

[4] Steve Ports an OFDM modem from Octave to C, https://www.rowetel.com/?p=5824

Quickstart
----------

Built as part of codec2-dev, see README for build instructions.

1. Generate 10 seconds of test frame bits, modulate, and play audio
   out of sound device (SoX v14.4.2):

     build_linux/src$ ./ofdm_mod --in /dev/zero --testframes 10 | play --type s16 --rate 8000 --channels 2 -

2. Generate 10 seconds of uncoded test frame bits, modulate, demodulate, count errors:

     build_linux/src$ ./ofdm_mod --in /dev/zero --testframes 10 | ./ofdm_demod --in /dev/null --testframes --log demod_dump.txt
     
   Use Octave to look at plots of C modem operation:

     $ cd ../../octave
     $ octave --no-gui
     octave:1> ofdm_demod_c("../build_linux/src/demod_dump.txt")

4. Run Octave versions of mod and demod (called tx and rx to avoid
   namespace clashes in Octave):

     $ cd ~/octave
     $ octave --no-gui
     octave:1> ofdm_tx("ofdm_test.raw","700D",10)
     octave:1> ofdm_rx("ofdm_test.raw")

   The Octave modulator ofdm_tx can simulate channel impairments, for
   example AWGN noise at an Eb/No of 4dB measured on 1400 bit/s raw
   uncoded bits:

     octave:1> ofdm_tx("ofdm_test.raw", "700D", 10, 4)

   The Octave versions use the same test frames as C so can interoperate.

     build_linux/src$ ./ofdm_demod --in ../../octave/ofdm_test.raw --out /dev/null --testframes --verbose 1

5. Run mod/demod with LDPC FEC; 4 frame interleaver, 60 seconds, 3dB
   Eb/No, Eb/No measured on 700 bit/s payload data bits.  For rate 1/2
   code this is equivalent to 0dB on 1400 bit/s uncoded bits (0dB
   Eb/No argument for ofdm_tx())

     octave:6> ofdm_ldpc_tx('ofdm_test.raw',"700D",4,60,3)
     octave:7> ofdm_ldpc_rx('ofdm_test.raw',"700D",4)

   C demodulator/LDPC decoder:
   
     build_linux/src$ ./ofdm_demod --in ../../octave/ofdm_test.raw --out /dev/null --verbose 1 --testframes --ldpc 1 --interleave 4

6. Run C mod/demod with LDPC and 2 frames interleaving:

     build_linux/src$ ./ofdm_mod --in /dev/zero --ldpc --testframes 60 --interleave 2 | ./ofdm_demod --out /dev/null --testframes --ldpc 1 --verbose 1 --interleave 2

7. Pass Codec 2 700C compressed speech through OFDM modem:

     build_linux/src$ ./c2enc 700C ../../raw/ve9qrp_10s.raw - --bitperchar | ./ofdm_mod --ldpc 1 --interleave 2 | ./ofdm_demod --ldpc 1 --interleave 2 | ./c2dec 700C - - --bitperchar | play --type s16 --rate 8000 --channels 1 -

8. Listen to signal through simulated fading channel in C:

   build_linux/src$ ./c2enc 700C ../../raw/ve9qrp_10s.raw - --bitperchar | ./ofdm_mod --ldpc 1 --interleave 4 | ./cohpsk_ch - - -20 -Fs 8000 --slow -f -5 | aplay -f S16

9. Run test frames through simulated channel in C:

   build_linux/src$ ./ofdm_mod --in /dev/zero --ldpc 1 --testframes 20 | ./cohpsk_ch - - -24 --Fs 8000 -f -10 --fast | ./ofdm_demod --out /dev/null --testframes --verbose 1 --ldpc 1

10. Run codec voice through simulated fast fading channel, just where it starts to fall over: 

   build_linux/src$ ./c2enc 700C ../../raw/ve9qrp.raw - --bitperchar | ./ofdm_mod --ldpc 1 --interleave 8 | ./cohpsk_ch - - -24 --Fs 8000 -f -10 --fast | ./ofdm_demod --ldpc 1 --verbose 1 --interleave 8 | ./c2dec 700C - - --bitperchar | aplay -f S16

11. FreeDV 1600 on the same channel conditions, roughly same quality at 8dB higher SNR:

   build_linux/src$ ./freedv_tx 1600 ../../raw/ve9qrp_10s.raw - - | ./cohpsk_ch - - -30 --Fs 8000 -f -10 --fast | ./freedv_rx 1600 - -  | aplay -f S16

12. Using FreeDV API test programs:

  build_linux/src$ ./freedv_tx 700D ../../raw/hts1a.raw - --testframes | ./freedv_rx 700D - /dev/null --testframes
  build_linux/src$ ./freedv_tx 700D ../../raw/hts1a.raw - | ./freedv_rx 700D - - | aplay -f S16

  With long interleaver times use a longer source file to allow interleaver time to sync, especially on poor
  channels:
  
  build_linux/src$ ./freedv_tx 700D ../../raw/ve9qrp.raw - - --interleave 8 | ./cohpsk_ch - - -26 --Fs 8000 -f -10 --fast | ./freedv_rx 700D - - --interleave 8 | aplay -f S16
  
FreeDV 2020 extensions
----------------------

13. 37 Carrier waveform with a (504,396) code:

   build_linux/src$ nc=37; ./ofdm_mod --in /dev/zero --testframes 300 --nc $nc --ldpc 2 --verbose 1 | ./cohpsk_ch - - -22.5 --Fs 8000 -f 10 --ssbfilt 1 | ./ofdm_demod --out /dev/null --testframes --nc $nc --verbose 1 --ldpc 2

   SNR3k(dB):  4.05 C/No: 38.8 PAPR: 10.8 
   BER......: 0.0348 Tbits: 1044792 Terrs: 36345
   Coded BER: 0.0094 Tbits: 820908 Terrs:  7717
   
14. 20.5ms symbol period, 31 carrier waveform, (504,396) code, but only
    312 data bits used, so we don't send unused data bits.  This means
    we need less carriers (so more power per carrier), and code rate
    is increased slightly (sorta).  Anyhoo, it works about 1.7dB
    better:

    build_linux/src$ nc=31; ./ofdm_mod --in /dev/zero --testframes 300 --ts 0.0205 --nc $nc --ldpc 2 --verbose 1 -p 312 | ./cohpsk_ch - - -21.6 --Fs 8000 -f 10 --ssbfilt 1 | ./ofdm_demod --out /dev/null --testframes --ts 0.0205 --nc $nc --verbose 1 --ldpc 2 -p 312

    SNR3k(dB):  2.21 C/No: 37.0 PAPR:  9.6 
    BER......: 0.0505 Tbits: 874020 Terrs: 44148
    Coded BER: 0.0096 Tbits: 649272 Terrs:  6230

15. Acquisition tests:

    Acquisition (getting sync) can be problematic in fading channels.
    Some special tests have been developed, that measure acquisition
    time on off air 700D samples at different time offsets:
    
    octave:61> ofdm_ldpc_rx("../wav/vk2tpm_004.wav", "700D", 1, "", 5, 4)
    build_linux/src$ ./ofdm_demod --in ../../wav/vk2tpm_004.wav --out /dev/null --verbose 2 --ldpc 1 --start_secs 5 --len_secs 4

    Different time offsets effectively tests the ability to sync on
    fading channel in different states.

    Stats for a series of these tests can be obtained with:
    
    octave:61> ofdm_time_sync("../wav/vk2tpm_004.wav", 30)
    <snip>pass: 30 fails: 0 mean: 1.35 var 0.51
    
Octave Acceptance Tests
-----------------------

Here are some useful tests for the Octave, uncoded modem.

The rate 1/2 LDPC code can correct up to about 10% raw BER, so a good
test is to run the modem at Eb/No operating points that produce just
less that BER=0.1

The BER2 measure truncates the effect of any start up transients,
e.g. as the frequency offset is tracked out.

1/ HF Multipath:

  octave:580> ofdm_tx("ofdm_test.raw","700D",60,4,'hf',20,-0.1)
  octave:581> ofdm_rx("ofdm_test.raw")
  BER2.: 0.0997 Tbits: 93752 Terrs:  9344

2/ AWGN:

  octave:582> ofdm_tx("ofdm_test.raw","700D",60,0,'awgn')
  octave:583> ofdm_rx("ofdm_test.raw")
  BER2.: 0.0827 Tbits: 96846 Terrs:  8008

C Acceptance Tests
------------------

Here are some useful tests for the LDPC coded C version of the modem, useful to verify any changes.

1/ AWGN channel, -2dB:

./ofdm_mod --in /dev/zero --ldpc 1 --testframes 60 --txbpf | ./cohpsk_ch - - -20 --Fs 8000 -f -10 | ./ofdm_demod --out /dev/null --testframes --verbose 1 --ldpc 1

SVN Rev 3671:
  SNR3k(dB): -1.85 C/No: 32.9 PAPR:  9.8
  BER......: 0.0815 Tbits: 98532 Terrs:  8031
  Coded BER: 0.0034 Tbits: 46368 Terrs:   157

2/ Fading HF channel:

./ofdm_mod --in /dev/zero --ldpc 1 --testframes 60 --txbpf | ./cohpsk_ch - - -24 --Fs 8000 -f -10 --fast | ./ofdm_demod --out /dev/null --testframes --verbose 1 --ldpc 1

SVN Rev 3671:
  SNR3k(dB):  2.15 C/No: 36.9 PAPR:  9.8
  BER......: 0.1015 Tbits: 88774 Terrs:  9012
  Coded BER: 0.0445 Tbits: 41776 Terrs:  1860

Note: 10% Raw BER operating point on both channels, as per design.  To
get a coded BER of around 1% on fading channels we need an interleaver
of about 2x the fade duration (e.g. --interleaver 16), above examples
are with interleaver == 1, so coded performance is not that great.

C Code
------

ofdm.c             - OFDM library
codec2_ofdm.h      - API header file for OFDM library 
ofdm_get_test_bits - generate OFDM test frames
ofdm_mod           - OFDM modulator command line program
ofdm_demod         - OFDM demodulator command line program, supports uncoded
                     (raw) and LDPC coded test frames, LDPC decoding
                     of codec data, interleaving, and can output LLRs to external
                     LDPC decoder                    
ofdm_put_test_bits - measure BER in OFDM test frames
unittest/tofdm     - Run C port of modem to compare with octave version (see octave/tofdm)
cohpsk_ch          - From COHPSK modem development, useful C channel simulator

Octave Scripts
--------------

ofdm_lib     - OFDM library 
ofdm_dev     - used for modem development, run various simulations
ofdm_tx      - modulate test frames to a file of sample, cam add channel impairments
ofdm_rx      - demod from a sample file and count errors
tofdm        - Compares Octave and C ports of modem
ofdm_ldpc_tx - OFDM modulator with interleaver and (224,112) LDPC code
ofdm_ldpc_rx - OFDM demodulator with interleaver and (224,112) LDPC code, interleaver sync

Specifications (Nominal FreeDV 700D configuration)
--------------------------------------------------

Modem.........: OFDM, pilot assisted coherent QPSK
Payload bits/s: 700
Text bits/s...: 25
Unique Word...: 10 bits
Carriers......: 17
RF bandwidth..: 944 Hz
Symbol period.: 18ms
Cyclic Prefix.: 2ms (note 1)
Pilot rate....: 1 in every 8 symbols
Frame Period..: 160ms
FEC...........: rate 1/2 (224,112) LDPC
Interleaving..: adjustable, suggest 1,4,8 or 16 frames
Operating point
  AWGN........: Eb/No -0.5dB SNR(3000Hz): -2.5dB (note 2)
  HF Multipath: Eb/No  4.0dB SNR(3000Hz):  2.0dB (note 3)
  
Freq offset.......: +/- 20  Hz   (sync range)
Freq drift........: +/- 0.2 Hz/s (for 0.5 dB loss)
Sample clock error: 1000 ppm

Notes:
  1/ Modem can cope with 2ms of multipath
  2/ Ideal SNR(3000) = Eb/No + 10*log10(Rb/B)
                     = -1 + 10*log10(1400/3000)
                     = -4.3 dB,
     So we have about 1.8dB overhead for synchronisation, implementation loss,
     and the text channel.
  3/ "CCIR Poor" HF Multipath channel used for testing is two path, 1Hz Doppler,
     1ms delay.
