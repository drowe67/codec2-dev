/*---------------------------------------------------------------------------*\

  FILE........: horus_l2.c
  AUTHOR......: David Rowe
  DATE CREATED: Dec 2015

  Horus telemetry layer 2 processing.  Takes an array of 8 bit payload
  data, generates parity bits for (23,12) Golay, interleaves data and
  parity bits, pre-pends a Unique Word for modem sync.  Caller is
  responsible for providing storage for output packet.

  [X] work out number of golay codewords rqd
  [X] function to return storgage rqd for output packet
  [X] code (not table based) based golay encoder
  [ ] unit test to check if working on payload
  [ ] code based interleaver
  [ ] unit test to run with/without errors
  [ ] conditional defines so just encoder on the target
  [ ] Octave code
      [ ] new protocol decoder
      [ ] golay in Octave
      [ ] test file based
      [ ] test real time

  Notes to Mark:
  [ ] dont send as RS232 (ie no start/stop bits)
  [ ] what are the word lengths of ints and longs on the uC?
  [ ] does it run fast enough on target?
  [ ] alternate legacy RTTY/new protocol
  [ ] test real time decode of both

  To test:

     src$ gcc golay23.c -o golay23 -Wall -DGOLAY23_UNITTEST -DRUN_TIME_TABLES
     src$ ./golay23

  To generate tables:
     src$ gcc golay23.c -o golay23 -Wall -DGOLAY23_MAKETABLES -DRUN_TIME_TABLES

\*---------------------------------------------------------------------------*/

#define RUN_TIME_TABLES

#include <assert.h>
#include <stdio.h>
#include <string.h>
#include "horus_l2.h"

static char uw[] = {'$','$'};

int get_syndrome(int pattern);
void golay23_init(void);
int golay23_decode(int received_codeword);

int horus_l2_get_num_tx_data_bytes(num_payload_data_bytes) {
    int num_payload_data_bits, num_golay_codewords;
    int num_tx_data_bits, num_tx_data_bytes;
    
    num_payload_data_bits = num_payload_data_bytes*8;
    num_golay_codewords = num_payload_data_bits/12;
    if (num_payload_data_bits % 12)
        num_golay_codewords++;

    num_tx_data_bits = sizeof(uw)*8 + num_payload_data_bits + num_golay_codewords*11;
    num_tx_data_bytes = num_tx_data_bits/8;
    if (num_tx_data_bits % 8)
        num_tx_data_bytes++;
    
    fprintf(stderr, "num_payload_data_bytes: %d\n", num_payload_data_bytes);
    fprintf(stderr, "num_golay_codewords...: %d\n", num_golay_codewords);
    fprintf(stderr, "num_tx_data_bits......: %d\n", num_tx_data_bits);
    fprintf(stderr, "num_tx_data_bytes.....: %d\n\n", num_tx_data_bytes);

    return num_tx_data_bytes;
}


/*
  The encoder will run on the payload on a small 8-bit uC.  As we are
  memory constrained so we do a lot of burrowing for bits out of
  packed arrays, and don't use a LUT for Golay encoding.  Hopefully it
  will run fast enough.
 */

int horus_l2_encode_tx_packet(unsigned char *output_tx_data,
                              unsigned char *input_payload_data,
                              int            num_payload_data_bytes)
{
    int            num_tx_data_bytes, num_payload_data_bits;
    unsigned char *pout = output_tx_data;
    int            ninbit, ingolay, ningolay, paritybyte, nparitybits;
    int            ninbyte, shift, inbit, golayparity, golayparitybit, i;

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
        fprintf(stderr, "inbit %d ninbyte: %d inbyte: 0x%02x inbit: %d\n", 
                ninbit, ninbyte, input_payload_data[ninbyte], inbit);
        ninbit++;

        /* build up input golay codeword */

        ingolay = ingolay | inbit;
        fprintf(stderr, "  ningolay: %d ingolay: 0x%04x\n", ningolay, ingolay);
        ningolay++;

        /* when we get 12 bits do a Golay encode */

        if (ningolay % 12) {
            ingolay <<= 1;
        }
        else {
            golayparity = get_syndrome(ingolay<<11);
            ingolay = 0;

            fprintf(stderr, "  golayparity: 0x%04x\n", golayparity);

            /* write parity bits to output data */

            for (i=0; i<11; i++) {
                golayparitybit = (golayparity >> (10-i)) & 0x1;
                paritybyte = paritybyte | golayparitybit;
                fprintf(stderr, "    i: %d golayparitybit: %d paritybyte: 0x%02x\n", 
                        i, golayparitybit, paritybyte);
                nparitybits++;
                if (nparitybits % 8) {
                   paritybyte <<= 1;
                }
                else {
                    /* OK we have a full byte ready */
                    *pout = paritybyte;
                    fprintf(stderr,"      Write paritybyte: 0x%02x\n", paritybyte);
                    pout++;
                    paritybyte = 0;
                }
            }
        }
    } /* while(.... */


    /* Complete final Golay encode, we may have partially finished ingolay, paritybyte */

    fprintf(stderr, "finishing up .....\n");

    if (ningolay % 12) {
        golayparity = get_syndrome(ingolay<<12);
        fprintf(stderr, "  golayparity: 0x%04x\n", golayparity);

        /* write parity bits to output data */

        for (i=0; i<11; i++) {
            golayparitybit = (golayparity >> (10 - i)) & 0x1;
            paritybyte = paritybyte | golayparitybit;
            fprintf(stderr, "    i: %d golayparitybit: %d paritybyte: 0x%02x\n", 
                    i, golayparitybit, paritybyte);
            nparitybits++;
            if (nparitybits % 8) {
                paritybyte <<= 1;
            }
            else {
                /* OK we have a full byte ready */
                *pout++ = (unsigned char)paritybyte;
                fprintf(stderr,"      Write paritybyte: 0x%02x\n", paritybyte);
                paritybyte = 0;
            }
        }
    }
 
    /* and final, partially complete, parity byte */

    if (nparitybits % 8) {
        paritybyte <<= 7 - (nparitybits % 8);  // use MS bits first
        *pout++ = (unsigned char)paritybyte;
        fprintf(stderr,"      Write last paritybyte: 0x%02x nparitybits: %d \n", paritybyte, nparitybits);
    }

    fprintf(stderr, "\npout - output_tx_data: %ld num_tx_data_bytes: %d\n",
            pout - output_tx_data, num_tx_data_bytes);
    assert(pout == (output_tx_data + num_tx_data_bytes));

    return num_tx_data_bytes;
}


void horus_l2_decode_rx_packet(unsigned char *output_payload_data,
                               unsigned char *input_rx_data,
                               int            num_payload_data_bytes)
{
    int            num_tx_data_bytes, num_payload_data_bits;
    unsigned char *pout = output_payload_data;
    unsigned char *pin  = input_rx_data;
    int            ninbit, ingolay, ningolay, paritybyte, nparitybits;
    int            ninbyte, shift, inbit, golayparitybit, i;

    num_tx_data_bytes = horus_l2_get_num_tx_data_bytes(num_payload_data_bytes);
    pin = input_rx_data + sizeof(uw) + num_payload_data_bytes;

    /* Read input data bits one at a time.  When we have 12 read 11 parity bits. Golay decode.
       Write decoded (output data) bits every time we have 8 of them. */

    num_payload_data_bits = num_payload_data_bytes*8;
    ninbit = 0;
    ingolay = 0;
    ningolay = 0;
    nparitybits = 0;
    paritybyte = *pin++;
    fprintf(stderr,"  Read paritybyte:0x%02x\n", paritybyte);
    pout = output_payload_data;

    while (ninbit < num_payload_data_bits) {

        /* extract input data bit */

        ninbyte = ninbit/8 + sizeof(uw);
        shift = 7 - (ninbit % 8);
        inbit = (input_rx_data[ninbyte] >> shift) & 0x1;
        fprintf(stderr, "inbit %d ninbyte: %d inbyte: 0x%02x inbit: %d\n", 
                ninbit, ninbyte, input_rx_data[ninbyte], inbit);
        ninbit++;

        /* build up golay codeword */

        ingolay = ingolay | inbit;
        fprintf(stderr, "  ningolay: %d ingolay: 0x%04x\n", ningolay, ingolay);
        ningolay++;
        ingolay <<= 1;

        /* when we get 12 data bits read parity bits */

        if ((ningolay % 12) == 0) {
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
                    fprintf(stderr,"  Read paritybyte: 0x%02x\n", paritybyte);
                }
            }

            fprintf(stderr, "  golay code word: 0x%04x\n", ingolay);
            fprintf(stderr, "  golay decode...: 0x%04x\n", golay23_decode(ingolay));
           
            /* Golay decode */

            
            /* write decoded/error corrected bits to output data */

            ingolay = 0;
        }
    } /* while(.... */


    fprintf(stderr, "finishing up .....\n");

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
                fprintf(stderr,"  Read paritybyte: 0x%02x\n", paritybyte);
            }
        }

        int codeword = (ingolay<<12) + golayparity;
        fprintf(stderr, "  golay code word: 0x%04x\n", codeword);
        fprintf(stderr, "  golay decode...: 0x%04x\n", golay23_decode(codeword));
    }
}


/* unit test designd to run on a Host PC */

int main(void) {
    unsigned char input_payload[] = {0x1,0x2,0x3};
    int num_tx_data_bytes = horus_l2_get_num_tx_data_bytes(sizeof(input_payload));
    unsigned char tx[num_tx_data_bytes];
    unsigned char output_payload[sizeof(input_payload)];
   int i;

    horus_l2_encode_tx_packet(tx, input_payload, sizeof(input_payload));

    fprintf(stderr, "\nTx Data:\n");
    for(i=0; i<num_tx_data_bytes; i++)
        fprintf(stderr, "  %02d 0x%02x\n", i, tx[i]);

    golay23_init();
    horus_l2_decode_rx_packet(output_payload, tx, sizeof(input_payload));

    return 0;
}


/*---------------------------------------------------------------------------*\

                                   GOLAY FUNCTIONS

\*---------------------------------------------------------------------------*/

/* File:    golay23.c
 * Title:   Encoder/decoder for a binary (23,12,7) Golay code
 * Author:  Robert Morelos-Zaragoza (robert@spectra.eng.hawaii.edu)
 * Date:    August 1994
 *
 * The binary (23,12,7) Golay code is an example of a perfect code, that is,
 * the number of syndromes equals the number of correctable error patterns.
 * The minimum distance is 7, so all error patterns of Hamming weight up to
 * 3 can be corrected. The total number of these error patterns is:
 *
 *       Number of errors         Number of patterns
 *       ----------------         ------------------
 *              0                         1
 *              1                        23
 *              2                       253
 *              3                      1771
 *                                     ----
 *    Total number of error patterns = 2048 = 2^{11} = number of syndromes
 *                                               --
 *                number of redundant bits -------^
 *
 * Because of its relatively low length (23), dimension (12) and number of
 * redundant bits (11), the binary (23,12,7) Golay code can be encoded and
 * decoded simply by using look-up tables. The program below uses a 16K
 * encoding table and an 8K decoding table.
 *
 * For more information, suggestions, or other ideas on implementing error
 * correcting codes, please contact me at (I'm temporarily in Japan, but
 * below is my U.S. address):
 *
 *                    Robert Morelos-Zaragoza
 *                    770 S. Post Oak Ln. #200
 *                      Houston, Texas 77056
 *
 *             email: robert@spectra.eng.hawaii.edu
 *
 *       Homework: Add an overall parity-check bit to get the (24,12,8)
 *                 extended Golay code.
 *
 * COPYRIGHT NOTICE: This computer program is free for non-commercial purposes.
 * You may implement this program for any non-commercial application. You may
 * also implement this program for commercial purposes, provided that you
 * obtain my written permission. Any modification of this program is covered
 * by this copyright.
 *
 * ==   Copyright (c) 1994  Robert Morelos-Zaragoza. All rights reserved.   ==
 */

#include "golay23.h"

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#define X22             0x00400000   /* vector representation of X^{22} */
#define X11             0x00000800   /* vector representation of X^{11} */
#define MASK12          0xfffff800   /* auxiliary vector for testing */
#define GENPOL          0x00000c75   /* generator polinomial, g(x) */

/* Global variables:
 *
 * pattern = error pattern, or information, or received vector
 * encoding_table[] = encoding table
 * decoding_table[] = decoding table
 * data = information bits, i(x)
 * codeword = code bits = x^{11}i(x) + (x^{11}i(x) mod g(x))
 * numerr = number of errors = Hamming weight of error polynomial e(x)
 * position[] = error positions in the vector representation of e(x)
 * recd = representation of corrupted received polynomial r(x) = c(x) + e(x)
 * decerror = number of decoding errors
 * a[] = auxiliary array to generate correctable error patterns
 */

static int inited =  0;

#ifdef RUN_TIME_TABLES
static int encoding_table[4096], decoding_table[2048];
#else
#include "golayenctable.h"
#include "golaydectable.h"
#endif

#ifdef GOLAY23_UNITTEST
static int position[23] = { 0x00000001, 0x00000002, 0x00000004, 0x00000008,
                            0x00000010, 0x00000020, 0x00000040, 0x00000080,
                            0x00000100, 0x00000200, 0x00000400, 0x00000800,
                            0x00001000, 0x00002000, 0x00004000, 0x00008000,
                            0x00010000, 0x00020000, 0x00040000, 0x00080000,
                            0x00100000, 0x00200000, 0x00400000 };
#endif

#ifdef RUN_TIME_TABLES
static int arr2int(int a[], int r)
/*
 * Convert a binary vector of Hamming weight r, and nonzero positions in
 * array a[1]...a[r], to a long integer \sum_{i=1}^r 2^{a[i]-1}.
 */
{
   int i;
   long mul, result = 0, temp;

   for (i=1; i<=r; i++) {
      mul = 1;
      temp = a[i]-1;
      while (temp--)
         mul = mul << 1;
      result += mul;
      }
   return(result);
}
#endif

void nextcomb(int n, int r, int a[])
/*
 * Calculate next r-combination of an n-set.
 */
{
  int  i, j;

  a[r]++;
  if (a[r] <= n)
      return;
  j = r - 1;
  while (a[j] == n - r + j)
     j--;
  for (i = r; i >= j; i--)
      a[i] = a[j] + i - j + 1;
  return;
}

int get_syndrome(int pattern)
/*
 * Compute the syndrome corresponding to the given pattern, i.e., the
 * remainder after dividing the pattern (when considering it as the vector
 * representation of a polynomial) by the generator polynomial, GENPOL.
 * In the program this pattern has several meanings: (1) pattern = infomation
 * bits, when constructing the encoding table; (2) pattern = error pattern,
 * when constructing the decoding table; and (3) pattern = received vector, to
 * obtain its syndrome in decoding.
 */
{
    int aux = X22;

    if (pattern >= X11)
       while (pattern & MASK12) {
           while (!(aux & pattern))
              aux = aux >> 1;
           pattern ^= (aux/X11) * GENPOL;
           }
    return(pattern);
}

/*---------------------------------------------------------------------------*\

  FUNCTION....: golay23_init()
  AUTHOR......: David Rowe
  DATE CREATED: 3 March 2013

  Call this once when you start your program to init the Golay tables.

\*---------------------------------------------------------------------------*/

void golay23_init(void) {
#ifdef RUN_TIME_TABLES
   int  i;
   long temp;
   int  a[4];
   int  pattern;

   /*
    * ---------------------------------------------------------------------
    *                  Generate ENCODING TABLE
    *
    * An entry to the table is an information vector, a 32-bit integer,
    * whose 12 least significant positions are the information bits. The
    * resulting value is a codeword in the (23,12,7) Golay code: A 32-bit
    * integer whose 23 least significant bits are coded bits: Of these, the
    * 12 most significant bits are information bits and the 11 least
    * significant bits are redundant bits (systematic encoding).
    * ---------------------------------------------------------------------
    */
    for (pattern = 0; pattern < 4096; pattern++) {
        temp = pattern << 11;          /* multiply information by X^{11} */
        encoding_table[pattern] = temp + get_syndrome(temp);/* add redundancy */
        }

   /*
    * ---------------------------------------------------------------------
    *                  Generate DECODING TABLE
    *
    * An entry to the decoding table is a syndrome and the resulting value
    * is the most likely error pattern. First an error pattern is generated.
    * Then its syndrome is calculated and used as a pointer to the table
    * where the error pattern value is stored.
    * ---------------------------------------------------------------------
    *
    * (1) Error patterns of WEIGHT 1 (SINGLE ERRORS)
    */
    decoding_table[0] = 0;
    decoding_table[1] = 1;
    temp = 1;
    for (i=2; i<= 23; i++) {
        temp *= 2;
        decoding_table[get_syndrome(temp)] = temp;
        }
   /*
    * (2) Error patterns of WEIGHT 2 (DOUBLE ERRORS)
    */
    a[1] = 1; a[2] = 2;
    temp = arr2int(a,2);
    decoding_table[get_syndrome(temp)] = temp;
    for (i=1; i<253; i++) {
        nextcomb(23,2,a);
        temp = arr2int(a,2);
        decoding_table[get_syndrome(temp)] = temp;
        }
   /*
    * (3) Error patterns of WEIGHT 3 (TRIPLE ERRORS)
    */
    a[1] = 1; a[2] = 2; a[3] = 3;
    temp = arr2int(a,3);
    decoding_table[get_syndrome(temp)] = temp;
    for (i=1; i<1771; i++) {
        nextcomb(23,3,a);
        temp = arr2int(a,3);
        decoding_table[get_syndrome(temp)] = temp;
    }
#endif
    inited = 1;
}

/*---------------------------------------------------------------------------*\

  FUNCTION....: golay23_encode()
  AUTHOR......: David Rowe
  DATE CREATED: 3 March 2013

  Given 12 bits of data retiurns a 23 bit codeword for transmission
  over the channel.

\*---------------------------------------------------------------------------*/

int golay23_encode(int data) {
    assert(inited);
    assert(data <= 0xfff);

    //printf("data: 0x%x\n", data);
    return encoding_table[data];
}

/*---------------------------------------------------------------------------*\

  FUNCTION....: golay23_decode()
  AUTHOR......: David Rowe
  DATE CREATED: 3 March 2013

  Given a 23 bit received codeword, returns the 12 bit corrected data.

\*---------------------------------------------------------------------------*/

int golay23_decode(int received_codeword) {
    assert(inited);
    assert((received_codeword < (1<<23)) && (received_codeword >= 0));

    //printf("syndrome: 0x%x\n", get_syndrome(received_codeword));
    return received_codeword ^= decoding_table[get_syndrome(received_codeword)];
}

int golay23_count_errors(int recd_codeword, int corrected_codeword)
{
    int errors = 0;
    int diff, i;

    diff = recd_codeword ^ corrected_codeword;
    for(i=0; i<23; i++) {
        if (diff & 0x1)
            errors++;
        diff >>= 1;
    }

    return errors;
}

#ifdef GOLAY23_UNITTEST

static int golay23_test(int error_pattern) {
    int data;
    int codeword;
    int recd;
    int pattern;
    int decerror;
    int i, tests;

    decerror = 0;
    tests = 0;

    for (data = 0; data<(1<<12); data++) {

        codeword = golay23_encode(data);
        recd = codeword ^ error_pattern;
        recd = golay23_decode(recd);
        pattern = (recd ^ codeword) >> 11;
        for (i=0; i<12; i++)
            if (pattern & position[i])
                decerror++;
        if (decerror) {
            printf("data: 0x%x codeword: 0x%x recd: 0x%x\n", data, codeword, recd);
            printf("there were %d decoding errors\n", decerror);
            exit(1);
        }
        tests++;
    }

    return tests;
}

int main(void)
{
   int i;
   int  tests;
   int a[4];
   int error_pattern;

   golay23_init();

   /* ---------------------------------------------------------------------
    *                        Generate DATA
    * ---------------------------------------------------------------------
    */

    /* Test all combinations of data and 1,2 or 3 errors */

    tests = 0;
    error_pattern = 1;
    for (i=0; i< 23; i++) {
        //printf("error_pattern: 0x%x\n", error_pattern);
        tests += golay23_test(error_pattern);
        error_pattern *= 2;
    }
    printf("%d 1 bit error tests performed OK!\n", tests);

    tests = 0;
    a[1] = 1; a[2] = 2;
    error_pattern = arr2int(a,2);
    tests += golay23_test(error_pattern);
    for (i=1; i<253; i++) {
        nextcomb(23,2,a);
        error_pattern = arr2int(a,2);
        //printf("error_pattern: 0x%x\n", error_pattern);
        tests += golay23_test(error_pattern);
    }
    printf("%d 2 bit error tests performed OK!\n", tests);

    tests = 0;
    a[1] = 1; a[2] = 2; a[3] = 3;
    error_pattern = arr2int(a,3);
    tests += golay23_test(error_pattern);
    for (i=1; i<1771; i++) {
        nextcomb(23,3,a);
        error_pattern = arr2int(a,3);
        //printf("error_pattern: 0x%x\n", error_pattern);
        tests += golay23_test(error_pattern);
    }
    printf("%d 3 bit error tests performed OK!\n", tests);

    return 0;
}
#endif

#ifdef GOLAY23_MAKETABLES
int main(int argc, char *argv[]) {
    FILE *f;
    int   i;

    golay23_init();

    f=fopen("golayenctable.h","wt");
    assert(f != NULL);

    fprintf(f,"/* Generated by golay23.c -DGOLAY23_MAKETABLE */\n\n");
    fprintf(f,"const int static encoding_table[]={\n");

    for (i=0; i<4095; i++)
        fprintf(f,"  0x%x,\n", encoding_table[i]);
    fprintf(f, "  0x%x\n};\n", encoding_table[i]);
    fclose(f);

    f=fopen("golaydectable.h","wt");
    assert(f != NULL);

    fprintf(f,"/* Generated by golay23.c -DGOLAY23_MAKETABLE */\n\n");
    fprintf(f,"const int static decoding_table[]={\n");

    for (i=0; i<2047; i++)
        fprintf(f,"  0x%x,\n", decoding_table[i]);
    fprintf(f, "  0x%x\n};\n", decoding_table[i]);
    fclose(f);

    return 0;
}

#endif


