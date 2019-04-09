# README_cohpsk.md

## Introduction

## Quickstart

1. BER test in AWGN channel with just less that 2% average bit error rate:

   ```
    $ ./cohpsk_get_test_bits - 5600 | ./cohpsk_mod - - --nd | ./cohpsk_ch - - -30 0 0 1 | ./cohpsk_demod - - --nd | ./cohpsk_put_test_bits -
    <snip>
    SNR3k(dB):  3.41 C/No: 38.2 PAPR:  8.1 
    BER: 0.017 Nbits: 5264 Nerrors: 92

  ```
  
## References

## C Code

## Octave Scripts


