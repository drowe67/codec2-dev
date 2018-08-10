/*---------------------------------------------------------------------------*\

  FILE........: horus_l2.c
  AUTHOR......: David Rowe
  DATE CREATED: Dec 2015

  Horus telemetry layer 2 processing.  Takes an array of 8 bit payload
  data, generates parity bits for a (23,12) Golay code, interleaves
  data and parity bits, pre-pends a Unique Word for modem sync.
  Caller is responsible for providing storage for output packet.

  1/ Unit test on a PC:

     $ gcc horus_l2.c -o horus_l2 -Wall -DHORUS_L2_UNITTEST
     $ ./horus_l2

     test 0: 22 bytes of payload data BER: 0.00 errors: 0
     test 0: 22 bytes of payload data BER: 0.01 errors: 0
     test 0: 22 bytes of payload data BER: 0.05 errors: 0
     test 0: 22 bytes of payload data BER: 0.10 errors: 7
     
     This indicates it's correcting all channel errors for 22 bytes of
     payload data, at bit error rate (BER) of 0, 0.01, 0.05.  It falls
     over at a BER of 0.10 which is expected.

  2/ To build with just the tx function, ie for linking with the payload
  firmware:

    $ gcc horus_l2.c golay23.c -c -Wall
    
  By default the RX side is #ifdef-ed out, leaving the minimal amount
  of code for tx.

  3/ Generate some tx_bits as input for testing with fsk_horus:
 
    $ gcc horus_l2.c golay23.c -o horus_l2 -Wall -DGEN_TX_BITS -DSCRAMBLER
    $ ./horus_l2
    $ more ../octave/horus_tx_bits_binary.txt
   
  4/ Unit testing interleaver:

    $ gcc horus_l2.c golay23.c -o horus_l2 -Wall -DINTERLEAVER -DTEST_INTERLEAVER -DSCRAMBLER

  5/ Compile for use as decoder called by fsk_horus.m and fsk_horus_stream.m:

    $ gcc horus_l2.c golay23.c -o horus_l2 -Wall -DDEC_RX_BITS -DHORUS_L2_RX

\*---------------------------------------------------------------------------*/

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include "horus_l2.h"
#include "golay23.h"

#ifdef HORUS_L2_UNITTEST
#define HORUS_L2_RX
#endif

static char uw[] = {'$','$'};

/* Function Prototypes ------------------------------------------------*/

#ifdef INTERLEAVER
static void interleave(unsigned char *inout, int nbytes, int dir);
#endif
#ifdef SCRAMBLER
static void scramble(unsigned char *inout, int nbytes);
#endif

/* Functions ----------------------------------------------------------*/

/*
   We are using a Golay (23,12) code which has a codeword 23 bits
   long.  The tx packet format is:

      | Unique Word | payload data bits | parity bits |

   This function works out how much storage the caller of
   horus_l2_encode_tx_packet() will need to store the tx packet
 */

int horus_l2_get_num_tx_data_bytes(int num_payload_data_bytes) {
    int num_payload_data_bits, num_golay_codewords;
    int num_tx_data_bits, num_tx_data_bytes;
    
    num_payload_data_bits = num_payload_data_bytes*8;
    num_golay_codewords = num_payload_data_bits/12;
    if (num_payload_data_bits % 12) /* round up to 12 bits, may mean some unused bits */
        num_golay_codewords++;

    num_tx_data_bits = sizeof(uw)*8 + num_payload_data_bits + num_golay_codewords*11;
    num_tx_data_bytes = num_tx_data_bits/8;
    if (num_tx_data_bits % 8) /* round up to nearest byte, may mean some unused bits */
        num_tx_data_bytes++;
    
    #ifdef DEBUG0
    fprintf(stderr, "\nnum_payload_data_bytes: %d\n", num_payload_data_bytes);
    fprintf(stderr, "num_golay_codewords...: %d\n", num_golay_codewords);
    fprintf(stderr, "num_tx_data_bits......: %d\n", num_tx_data_bits);
    fprintf(stderr, "num_tx_data_bytes.....: %d\n\n", num_tx_data_bytes);
    #endif

    return num_tx_data_bytes;
}

void horus_l2_init(void) {
    golay23_init();
}

/*
  Takes an array of payload data bytes, prepends a unique word and appends
  parity bits.

  The encoder will run on the payload on a small 8-bit uC.  As we are
  memory constrained so we do a lot of burrowing for bits out of
  packed arrays, and don't use a LUT for Golay encoding.  Hopefully it
  will run fast enough.  This was quite difficult to get going,
  suspect there is a better way to write this.  Oh well, have to start
  somewhere.
 */

int horus_l2_encode_tx_packet(unsigned char *output_tx_data,
                              unsigned char *input_payload_data,
                              int            num_payload_data_bytes)
{
    int            num_tx_data_bytes, num_payload_data_bits;
    unsigned char *pout = output_tx_data;
    int            ninbit, ningolay, nparitybits;
    int32_t        ingolay, paritybyte, inbit, golayparity;
    int            ninbyte, shift, golayparitybit, i;

    num_tx_data_bytes = horus_l2_get_num_tx_data_bytes(num_payload_data_bytes);
    memcpy(pout, uw, sizeof(uw)); pout += sizeof(uw);
    memcpy(pout, input_payload_data, num_payload_data_bytes); pout += num_payload_data_bytes;

    /* Read input bits one at a time.  Fill input Golay codeword.  Find output Golay codeword.
       Write this to parity bits.  Write parity bytes when we have 8 parity bits.  Bits are
       written MSB first. */

    num_payload_data_bits = num_payload_data_bytes*8;
    ninbit = 0;
    ingolay = 0;
    ningolay = 0;
    paritybyte = 0;
    nparitybits = 0;

    while (ninbit < num_payload_data_bits) {

        /* extract input data bit */

        ninbyte = ninbit/8;
        shift = 7 - (ninbit % 8);
        inbit = (input_payload_data[ninbyte] >> shift) & 0x1;
        #ifdef DEBUG1
        fprintf(stderr, "inbit %d ninbyte: %d inbyte: 0x%02x inbit: %d\n", 
                ninbit, ninbyte, input_payload_data[ninbyte], inbit);
        #endif
        ninbit++;

        /* build up input golay codeword */

        ingolay = ingolay | inbit;
        ningolay++;

        /* when we get 12 bits do a Golay encode */

        if (ningolay % 12) {
            ingolay <<= 1;
        }
        else {
            #ifdef DEBUG0
            fprintf(stderr, "  ningolay: %d ingolay: 0x%04x\n", ningolay, ingolay);
            #endif
            golayparity = golay23_syndrome(ingolay<<11);
            ingolay = 0;

            #ifdef DEBUG0
            fprintf(stderr, "  golayparity: 0x%04x\n", golayparity);
            #endif

            /* write parity bits to output data */

            for (i=0; i<11; i++) {
                golayparitybit = (golayparity >> (10-i)) & 0x1;
                paritybyte = paritybyte | golayparitybit;
                #ifdef DEBUG0
                fprintf(stderr, "    i: %d golayparitybit: %d paritybyte: 0x%02x\n", 
                        i, golayparitybit, paritybyte);
                #endif
                nparitybits++;
                if (nparitybits % 8) {
                   paritybyte <<= 1;
                }
                else {
                    /* OK we have a full byte ready */
                    *pout = paritybyte;
                    #ifdef DEBUG0
                    fprintf(stderr,"      Write paritybyte: 0x%02x\n", paritybyte);
                    #endif
                    pout++;
                    paritybyte = 0;
                }
            }
        }
    } /* while(.... */


    /* Complete final Golay encode, we may have partially finished ingolay, paritybyte */

    #ifdef DEBUG0
    fprintf(stderr, "finishing up .....\n");
    #endif

    if (ningolay % 12) {
        ingolay >>= 1;
        golayparity = golay23_syndrome(ingolay<<12);
        #ifdef DEBUG0
        fprintf(stderr, "  ningolay: %d ingolay: 0x%04x\n", ningolay, ingolay);
        fprintf(stderr, "  golayparity: 0x%04x\n", golayparity);
        #endif

        /* write parity bits to output data */

        for (i=0; i<11; i++) {
            golayparitybit = (golayparity >> (10 - i)) & 0x1;
            paritybyte = paritybyte | golayparitybit;
            #ifdef DEBUG1
            fprintf(stderr, "    i: %d golayparitybit: %d paritybyte: 0x%02x\n", 
                    i, golayparitybit, paritybyte);
            #endif
            nparitybits++;
            if (nparitybits % 8) {
                paritybyte <<= 1;
            }
            else {
                /* OK we have a full byte ready */
                *pout++ = (unsigned char)paritybyte;
                #ifdef DEBUG0
                fprintf(stderr,"      Write paritybyte: 0x%02x\n", paritybyte);
                #endif
                paritybyte = 0;
            }
        }
    }
 
    /* and final, partially complete, parity byte */

    if (nparitybits % 8) {
        paritybyte <<= 7 - (nparitybits % 8);  // use MS bits first
        *pout++ = (unsigned char)paritybyte;
        #ifdef DEBUG0
        fprintf(stderr,"      Write last paritybyte: 0x%02x nparitybits: %d \n", paritybyte, nparitybits);
        #endif
    }

    #ifdef DEBUG0
    fprintf(stderr, "\npout - output_tx_data: %ld num_tx_data_bytes: %d\n",
            pout - output_tx_data, num_tx_data_bytes);
    #endif
    assert(pout == (output_tx_data + num_tx_data_bytes));

    /* optional interleaver - we dont interleave UW */

    #ifdef INTERLEAVER
    interleave(&output_tx_data[sizeof(uw)], num_tx_data_bytes-2, 0);
    #endif

    /* optional scrambler to prevent long strings of the same symbol
       which upsets the modem - we dont scramble UW */

    #ifdef SCRAMBLER
    scramble(&output_tx_data[sizeof(uw)], num_tx_data_bytes-2);
    #endif

    return num_tx_data_bytes;
}


#ifdef HORUS_L2_RX
void horus_l2_decode_rx_packet(unsigned char *output_payload_data,
                               unsigned char *input_rx_data,
                               int            num_payload_data_bytes)
{
    int            num_payload_data_bits;
    unsigned char *pout = output_payload_data;
    unsigned char *pin  = input_rx_data;
    int            ninbit, ingolay, ningolay, paritybyte, nparitybits;
    int            ninbyte, shift, inbit, golayparitybit, i, outbit, outbyte, noutbits, outdata;
    #if defined(SCRAMBLER) || defined(INTERLEAVER)
    int num_tx_data_bytes = horus_l2_get_num_tx_data_bytes(num_payload_data_bytes);
    #endif
    
    /* optional scrambler and interleaver - we dont interleave UW */

    #ifdef SCRAMBLER
    scramble(&input_rx_data[sizeof(uw)], num_tx_data_bytes-2);
    #endif

    #ifdef INTERLEAVER
    interleave(&input_rx_data[sizeof(uw)], num_tx_data_bytes-2, 1);
    #endif

    pin = input_rx_data + sizeof(uw) + num_payload_data_bytes;

    /* Read input data bits one at a time.  When we have 12 read 11 parity bits. Golay decode.
       Write decoded (output data) bits every time we have 8 of them. */

    num_payload_data_bits = num_payload_data_bytes*8;
    ninbit = 0;
    ingolay = 0;
    ningolay = 0;
    nparitybits = 0;
    paritybyte = *pin++;
    #ifdef DEBUG0
    fprintf(stderr,"  Read paritybyte: 0x%02x\n", paritybyte);
    #endif
    pout = output_payload_data;
    noutbits = 0;
    outbyte = 0;

    while (ninbit < num_payload_data_bits) {

        /* extract input data bit */

        ninbyte = ninbit/8 + sizeof(uw);
        shift = 7 - (ninbit % 8);
        inbit = (input_rx_data[ninbyte] >> shift) & 0x1;
        #ifdef DEBUG1
        fprintf(stderr, "inbit %d ninbyte: %d inbyte: 0x%02x inbit: %d\n", 
                ninbit, ninbyte, input_rx_data[ninbyte], inbit);
        #endif
        ninbit++;

        /* build up golay codeword */

        ingolay = ingolay | inbit;
        ningolay++;
        ingolay <<= 1;

        /* when we get 12 data bits start reading parity bits */

        if ((ningolay % 12) == 0) {
            #ifdef DEBUG0
            fprintf(stderr, "  ningolay: %d ingolay: 0x%04x\n", ningolay, ingolay>>1);
            #endif
            for (i=0; i<11; i++) {
                shift = 7 - (nparitybits % 8);
                golayparitybit = (paritybyte >> shift) & 0x1;
                ingolay |= golayparitybit;
                if (i != 10)
                    ingolay <<=1;
                nparitybits++;
                if ((nparitybits % 8) == 0) {
                    /* OK grab a new byte */
                    paritybyte = *pin++;
                    #ifdef DEBUG0
                    fprintf(stderr,"  Read paritybyte: 0x%02x\n", paritybyte);
                    #endif
                }
            }

            #ifdef DEBUG0
            fprintf(stderr, "  golay code word: 0x%04x\n", ingolay);
            fprintf(stderr, "  golay decode...: 0x%04x\n", golay23_decode(ingolay));
            #endif
           
            /* write decoded/error corrected bits to output payload data */

            outdata = golay23_decode(ingolay) >> 11;
            #ifdef DEBUG0
            fprintf(stderr, "  outdata...: 0x%04x\n", outdata);
            #endif

            for(i=0; i<12; i++) {   
                shift = 11 - i;
                outbit = (outdata >> shift) & 0x1;
                outbyte |= outbit;
                noutbits++;
                if (noutbits % 8) {
                    outbyte <<= 1;
                }
                else {
                    #ifdef DEBUG0
                    fprintf(stderr, "  output payload byte: 0x%02x\n", outbyte);
                    #endif
                    *pout++ = outbyte;
                    outbyte = 0;
                }
            }

            ingolay = 0;
        }
    } /* while(.... */


    #ifdef DEBUG0
    fprintf(stderr, "finishing up .....\n");
    #endif

    /* Complete final Golay decode  */

    int golayparity = 0;
    if (ningolay % 12) {
        for (i=0; i<11; i++) {
            shift = 7 - (nparitybits % 8);
            golayparitybit = (paritybyte >> shift) & 0x1;
            golayparity |= golayparitybit;
            if (i != 10)
                golayparity <<=1;
            nparitybits++;
            if ((nparitybits % 8) == 0) {
                /* OK grab a new byte */
                paritybyte = *pin++;
                #ifdef DEBUG0
                fprintf(stderr,"  Read paritybyte: 0x%02x\n", paritybyte);
                #endif
            }
        }

        ingolay >>= 1;
        int codeword = (ingolay<<12) + golayparity;
        #ifdef DEBUG0
        fprintf(stderr, "  ningolay: %d ingolay: 0x%04x\n", ningolay, ingolay);
        fprintf(stderr, "  golay code word: 0x%04x\n", codeword);
        fprintf(stderr, "  golay decode...: 0x%04x\n", golay23_decode(codeword));
        #endif

        outdata = golay23_decode(codeword) >> 11;
        #ifdef DEBUG0
        fprintf(stderr, "  outdata...: 0x%04x\n", outdata);
        fprintf(stderr, "  num_payload_data_bits: %d noutbits: %d\n", num_payload_data_bits, noutbits);
        #endif

        /* write final byte */

        int ntogo = num_payload_data_bits - noutbits;
        for(i=0; i<ntogo; i++) {   
            shift = ntogo - i;
            outbit = (outdata >> shift) & 0x1;
            outbyte |= outbit;
            noutbits++;
            if (noutbits % 8) {
                outbyte <<= 1;
            }
            else {
                #ifdef DEBUG0
                fprintf(stderr, "  output payload byte: 0x%02x\n", outbyte);
                #endif
                *pout++ = outbyte;
                outbyte = 0;
            }
        }
    }

    #ifdef DEBUG0
    fprintf(stderr, "\npin - output_payload_data: %ld num_payload_data_bytes: %d\n",
            pout - output_payload_data, num_payload_data_bytes);
    #endif

    assert(pout == (output_payload_data + num_payload_data_bytes));

}
#endif

#ifdef INTERLEAVER

uint16_t primes[] = {
    2,      3,      5,      7,      11,     13,     17,     19,     23,     29, 
    31,     37,     41,     43,     47,     53,     59,     61,     67,     71, 
    73,     79,     83,     89,     97,     101,    103,    107,    109,    113, 
    127,    131,    137,    139,    149,    151,    157,    163,    167,    173, 
    179,    181,    191,    193,    197,    199,    211,    223,    227,    229, 
    233,    239,    241,    251,    257,    263,    269,    271,    277,    281, 
    283,    293,    307,    311,    313,    317,    331,    337,    347
};

void interleave(unsigned char *inout, int nbytes, int dir)
{
    /* note: to work on small uCs (e.g. AVR) needed to declare specific words sizes */
    uint16_t nbits = (uint16_t)nbytes*8;
    uint32_t i, j, n, ibit, ibyte, ishift, jbyte, jshift;
    uint32_t b;
    unsigned char out[nbytes];

    memset(out, 0, nbytes);
           
    /* b chosen to be co-prime with nbits, I'm cheating by just finding the 
       nearest prime to nbits.  It also uses storage, is run on every call,
       and has an upper limit.  Oh Well, still seems to interleave OK. */
    i = 1;
    uint16_t imax = sizeof(primes)/sizeof(uint16_t);
    while ((primes[i] < nbits) && (i < imax))
        i++;
    b = primes[i-1];

    for(n=0; n<nbits; n++) {

        /*
          "On the Analysis and Design of Good Algebraic Interleavers", Xie et al,eq (5)
        */

        i = n;
        j = (b*i) % nbits; /* note these all need to be 32-bit ints to make multiply work without overflow */
        
        if (dir) {
            uint16_t tmp = j;
            j = i;
            i = tmp;
        }
        
        #ifdef DEBUG0
        printf("i: %d j: %d\n",i, j);
        #endif

        /* read bit i and write to bit j postion */

        ibyte = i/8;
        ishift = i%8;
        ibit = (inout[ibyte] >> ishift) & 0x1;

        jbyte = j/8;
        jshift = j%8;

        /* write jbit to ibit position */

        out[jbyte] |= ibit << jshift; // replace with i-th bit
        //out[ibyte] |= ibit << ishift; // replace with i-th bit
    }
 
    memcpy(inout, out, nbytes);

    #ifdef DEBUG0
    printf("\nInterleaver Out:\n");
    for (i=0; i<nbytes; i++)
        printf("%02d 0x%02x\n", i, inout[i]);
    #endif
}
#endif


#ifdef TEST_INTERLEAVER
int main(void) {
    int nbytes = 43;
    unsigned char inout[nbytes];
    unsigned char inter[nbytes];
    unsigned char incopy[nbytes];
    int i;

    /* copy of input for later comp   */

    for(i=0; i<nbytes; i++)
        inout[i] = incopy[i] = rand() & 0xff;    
    
    interleave(inout, nbytes, 0);    /* interleave                     */
    memcpy(inter, inout, nbytes);    /* snap shot of interleaved bytes */
    interleave(inout, nbytes, 1);    /* de-interleave                  */

    /* all ones in last col means it worked! */

    for(i=0; i<nbytes; i++) {
        printf("%d 0x%02x 0x%02x 0x%02x %d\n", 
               i, incopy[i], inter[i], inout[i],  incopy[i] == inout[i]);
        assert(incopy[i] == inout[i]);
    }
    printf("Interleaver tested OK!\n");

    return 0;
}
#endif


#ifdef SCRAMBLER

/* 16 bit DVB additive scrambler as per Wikpedia example */

void scramble(unsigned char *inout, int nbytes)
{
    int nbits = nbytes*8;
    int i, ibit, ibits, ibyte, ishift, mask;
    uint16_t scrambler = 0x4a80;  /* init additive scrambler at start of every frame */
    uint16_t scrambler_out;

    /* in place modification of each bit */

    for(i=0; i<nbits; i++) {

        scrambler_out = ((scrambler & 0x2) >> 1) ^ (scrambler & 0x1);

        /* modify i-th bit by xor-ing with scrambler output sequence */

        ibyte = i/8;
        ishift = i%8;
        ibit = (inout[ibyte] >> ishift) & 0x1;
        ibits = ibit ^ scrambler_out;                  // xor ibit with scrambler output

        mask = 1 << ishift;
        inout[ibyte] &= ~mask;                  // clear i-th bit
        inout[ibyte] |= ibits << ishift;         // set to scrambled value

        /* update scrambler */

        scrambler >>= 1;
        scrambler |= scrambler_out << 14;

        #ifdef DEBUG0
        printf("i: %02d ibyte: %d ishift: %d ibit: %d ibits: %d scrambler_out: %d\n", 
               i, ibyte, ishift, ibit, ibits, scrambler_out);
        #endif

    }

    #ifdef DEBUG0
    printf("\nScrambler Out:\n");
    for (i=0; i<nbytes; i++)
        printf("%02d 0x%02x\n", i, inout[i]);
    #endif
}
#endif

#ifdef HORUS_L2_UNITTEST

/*
  Test function to construct a packet of payload data, encode, add
  some bit errors, decode, count errors.
*/

int test_sending_bytes(int nbytes, float ber, int error_pattern) {
    unsigned char input_payload[nbytes];
    int num_tx_data_bytes = horus_l2_get_num_tx_data_bytes(sizeof(input_payload));
    unsigned char tx[num_tx_data_bytes];
    unsigned char output_payload[sizeof(input_payload)];
    int b, nbiterrors = 0;
    int i;

    for(i=0; i<nbytes; i++)
        input_payload[i] = i;

    horus_l2_encode_tx_packet(tx, input_payload, sizeof(input_payload));

    #ifdef DEBUG0
    fprintf(stderr, "\nTx Data:\n");
    for(i=0; i<num_tx_data_bytes; i++)
        fprintf(stderr, "  %02d 0x%02x\n", i, tx[i]);
    #endif

    /* insert random bit errors */

    if (error_pattern == 0) {
        float r;
        for(i=0; i<num_tx_data_bytes; i++) {
            for (b=0; b<8; b++) {
                r = (float)rand()/RAND_MAX;
                if (r < ber) {
                    unsigned char mask = (1<<b);
                    #ifdef DEBUG1
                    fprintf("mask: 0x%x tx[%d] = 0x%x ", mask, i, tx[i]);
                    #endif
                    tx[i] ^= mask;
                    #ifdef DEBUG1
                    fprintf("0x%x\n", tx[i]);
                    #endif
                    nbiterrors++;
                }
            }
        }
    }

    /* insert and error burst */

    if (error_pattern == 1) {
        tx[2] ^= 0xff;
        tx[3] ^= 0xff;
    }

    /* insert 1 error every 12 bits, this gives up to 3 errors per 23
       bit codeword, which is the limit of the code */

    if (error_pattern == 2) {
        int bn = 0;
        for(i=0; i<num_tx_data_bytes; i++) {
            for (b=0; b<8; b++) {
                bn++;
                if ((bn % 12) == 0) {
                    unsigned char mask = (1<<(7-b));
                    #ifdef DEBUG1
                    fprintf("mask: 0x%x tx[%d] = 0x%x ", mask, i, tx[i]);
                    #endif
                    tx[i] ^= mask;
                    #ifdef DEBUG1
                    fprintf("0x%x\n", tx[i]);
                    #endif
                    nbiterrors++;
                }
            }
        }
    }

    #ifdef DEBUG0
    fprintf(stderr, "\nTx Data after errors:\n");
    for(i=0; i<num_tx_data_bytes; i++)
        fprintf(stderr, "  %02d 0x%02x\n", i, tx[i]);
    #endif

    #ifdef DEBUG0
    fprintf(stderr, "nbiterrors: %d BER: %3.2f\n", nbiterrors, (float)nbiterrors/(num_tx_data_bytes*8));
    #endif

    golay23_init();
    horus_l2_decode_rx_packet(output_payload, tx, sizeof(input_payload));

    #ifdef DEBUG0
    fprintf(stderr, "\nOutput Payload:\n");
    for(i=0; i<sizeof(input_payload); i++)
        fprintf(stderr, "  %02d 0x%02x\n", i, output_payload[i]);
    #endif

    /* count bit errors */

    int nerr = 0;
    for(i=0; i<nbytes; i++) {
        int error_pattern = input_payload[i] ^ output_payload[i];
        for(b=0; b<8; b++)
            nerr += (error_pattern>>b) & 0x1;
    }
    
    return nerr;
}

/* unit test designed to run on a PC */

int main(void) {
    printf("test 0: BER: 0.00 ...........: %d\n", test_sending_bytes(22, 0.00, 0));
    printf("test 1: BER: 0.01 ...........: %d\n", test_sending_bytes(22, 0.01, 0));
    printf("test 2: BER: 0.05 ...........: %d\n", test_sending_bytes(22, 0.05, 0));

    /* we expect this always to fail, as chance of > 3 errors/codeword is high */

    printf("test 3: BER: 0.10 ...........: %d\n", test_sending_bytes(22, 0.10, 0));

    /* -DINTERLEAVER will make this puppy pass */

    printf("test 4: 8 bit burst error....: %d\n", test_sending_bytes(22, 0.00, 1));

    /* Insert 2 errors in every codeword, the maximum correction
       capability of a Golay (23,12) code. note this one will fail
       with -DINTERLEAVER, as we can't guarantee <= 3 errors per
       codeword after interleaving */

    printf("test 5: 1 error every 12 bits: %d\n", test_sending_bytes(22, 0.00, 2));
    return 0;
}
#endif

/* Horus binary packet */

struct TBinaryPacket
{
    uint8_t     PayloadID;
    uint16_t	Counter;
    uint8_t	Hours;
    uint8_t	Minutes;
    uint8_t	Seconds;
    float	Latitude;
    float	Longitude;
    uint16_t  	Altitude;
    uint8_t     Speed;       // Speed in Knots (1-255 knots)
    uint8_t     Sats;
    int8_t      Temp;        // Twos Complement Temp value.
    uint8_t     BattVoltage; // 0 = 0.5v, 255 = 2.0V, linear steps in-between.
    uint16_t    Checksum;    // CRC16-CCITT Checksum.
}  __attribute__ ((packed));

#ifdef GEN_TX_BITS
/* generate a file of tx_bits to modulate using fsk_horus.m for modem simulations */

int main(void) {
    int nbytes = sizeof(struct TBinaryPacket);
    struct TBinaryPacket input_payload;
    int num_tx_data_bytes = horus_l2_get_num_tx_data_bytes(nbytes);
    unsigned char tx[num_tx_data_bytes];
    int i;

    /* all zeros is nastiest sequence for demod before scrambling */

    memset(&input_payload, 0, nbytes);
    input_payload.Checksum = horus_l2_gen_crc16((unsigned char*)&input_payload, nbytes-2);

    horus_l2_encode_tx_packet(tx, (unsigned char*)&input_payload, nbytes);
    
    FILE *f = fopen("../octave/horus_tx_bits_binary.txt","wt");
    assert(f != NULL);
    int b, tx_bit;
    for(i=0; i<num_tx_data_bytes; i++) {
        for(b=0; b<8; b++) {
            tx_bit = (tx[i] >> (7-b)) & 0x1; /* msb first */
            fprintf(f,"%d ", tx_bit);
        }
    }
    fclose(f);

    return 0;
}
#endif


#ifdef DEC_RX_BITS

/* Decode a binary file rx_bytes, e.g. from fsk_horus.m */

int main(void) {
    int nbytes = 22;
    unsigned char output_payload[nbytes];
    int num_tx_data_bytes = horus_l2_get_num_tx_data_bytes(nbytes);

    /* real world data horus payload generated when running tx above */
    unsigned char rx[45] = {
        0x24,0x24,0x01,0x0b,0x00,0x00,0x05,0x3b,0xf2,0xa7,0x0b,0xc2,0x1b,
        0xaa,0x0a,0x43,0x7e,0x00,0x05,0x00,0x25,0xc0,0xce,0xbb,0x36,0x69,
        0x50,0x00,0x41,0xb0,0xa6,0x5e,0x91,0xa2,0xa3,0xf8,0x1d,0x00,0x00,
        0x0c,0x76,0xc6,0x05,0xb0,0xb8};
    int i, ret;

    assert(num_tx_data_bytes == 45);

    #define READ_FILE /* overwrite tx[] above, that's OK */
    #ifdef READ_FILE
    FILE *f = fopen("../octave/horus_rx_bits_binary.bin","rb");
    assert(f != NULL);
    ret = fread(rx, sizeof(char), num_tx_data_bytes, f);
    assert(ret == num_tx_data_bytes);
    fclose(f);
    #endif

    golay23_init();
    horus_l2_decode_rx_packet(output_payload, rx, nbytes);

    #ifdef HEX_DUMP
    fprintf(stderr, "\nOutput Payload:\n");
    for(i=0; i<nbytes; i++)
        fprintf(stderr, "  %02d 0x%02x 0x%02x\n", i, output_payload[i], rx[i+2]);
    #endif

    struct TBinaryPacket h;
    assert(sizeof(h) == nbytes);
    memcpy(&h, output_payload, nbytes);

    uint16_t crc_rx = horus_l2_gen_crc16(output_payload, nbytes-2);
    char crc_str[80];
    
    if (crc_rx == h.Checksum) {
        sprintf(crc_str, "CRC OK");
    } else {
        sprintf(crc_str, "CRC BAD");
    }

    fprintf(stderr, "%d,%d,%02d:%02d:%02d,%f,%f,%d,%d,%d,%d,%d,%04x %s\n",
        h.PayloadID, h.Counter, h.Hours, h.Minutes, h.Seconds,
        h.Latitude, h.Longitude, h.Altitude, h.Speed, h.Sats, h.Temp, 
            h.BattVoltage, h.Checksum, crc_str);
    
    /* Hex ASCII file output */

    #define WRITE_HEX_FILE /* overwrite tx[] above, that's OK */
    #ifdef WRITE_HEX_FILE
    FILE *fh = fopen("../octave/horus_rx_bits_hex.txt","wt");
    assert(fh != NULL);
    for(i=0; i<nbytes; i++) {
        fprintf(fh, "%02X", (unsigned int)output_payload[i]);
    }
    fclose(fh);
    #endif

    return 0;
}
#endif

// from http://stackoverflow.com/questions/10564491/function-to-calculate-a-crc16-checksum

unsigned short horus_l2_gen_crc16(unsigned char* data_p, unsigned char length) {
    unsigned char x;
    unsigned short crc = 0xFFFF;

    while (length--){
        x = crc >> 8 ^ *data_p++;
        x ^= x>>4;
        crc = (crc << 8) ^ ((unsigned short)(x << 12)) ^ ((unsigned short)(x <<5)) ^ ((unsigned short)x);
    }
    return crc;
}

