/*---------------------------------------------------------------------------*\

  FILE........: drs232_ldpc.c
  AUTHOR......: David Rowe
  DATE CREATED: Sep 2016

  Looks for a unique word in series of soft decision symbols.  When
  found, deframes a RS232 encoded frame of soft decision bit, LDPC
  decodes, and outputs a frame of packed bytes.  Used for high biit
  rate Horus SSTV reception.

  Frame format:

    16 bytes 0x55 - 0xabcdef01 UW - 256 bytes of payload - 2 bytes CRC - 65 bytes LPDC Parity

  Each byte is encoded as a 10 bit RS232 serial word: 
  
    0 LSB .... MSB 1

  Building:
   
    $ gcc drs232_ldpc.c mpdecode_core.c -o drs232_ldpc -Wall -lm

\*---------------------------------------------------------------------------*/


/*
  Copyright (C) 2016 David Rowe

  All rights reserved.

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU Lesser General Public License version 2.1, as
  published by the Free Software Foundation.  This program is
  distributed in the hope that it will be useful, but WITHOUT ANY
  WARRANTY; without even the implied warranty of MERCHANTABILITY or
  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public
  License for more details.

  You should have received a copy of the GNU Lesser General Public License
  along with this program; if not, see <http://www.gnu.org/licenses/>.
*/

#include <assert.h>
#include <stdlib.h>
#include <errno.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>

#include "mpdecode_core.h"

/* Machine generated consts, H_rows, H_cols, test input/output data to
   change LDPC code regenerate this file. */

#include "H2064_516_sparse.h"  

/* states -----------------------------------------------*/

#define LOOK_FOR_UW    0
#define COLLECT_PACKET 1

/* packet parameters */

#define UW_BYTES               4
#define UW_BITS                40
#define UW_ALLOWED_ERRORS      5
#define BYTES_PER_PACKET       256
#define CRC_BYTES              2
#define PARITY_BYTES           65
#define BITS_PER_BYTE          10
#define UNPACKED_PACKET_BYTES  ((UW_BYTES+BYTES_PER_PACKET+CRC_BYTES)*BITS_PER_BYTE)
#define SYMBOLS_PER_PACKET     (BYTES_PER_PACKET+CRC_BYTES+PARITY_BYTES)*BITS_PER_BYTE

/* UW pattern we look for, including start/stop bits */

uint8_t uw[] = {
    /* 0xb                0xa */
    0, 1, 1, 0, 1, 0, 1, 0, 1, 1,
    /* 0xd                0xc */
    0, 1, 0, 1, 1, 0, 0, 1, 1, 1,
    /* 0xf                0xe */
    0, 1, 1, 1, 1, 0, 1, 1, 1, 1,
    /* 0x1                0x0 */
    0, 1, 0, 0, 0, 0, 0, 0, 0, 1,
};


// from http://stackoverflow.com/questions/10564491/function-to-calculate-a-crc16-checksum

unsigned short gen_crc16(unsigned char* data_p, int length){
    unsigned char x;
    unsigned short crc = 0xFFFF;

    while (length--){
        x = crc >> 8 ^ *data_p++;
        x ^= x>>4;
        crc = (crc << 8) ^ ((unsigned short)(x << 12)) ^ ((unsigned short)(x <<5)) ^ ((unsigned short)x);
    }
    
    return crc;
}


int main(int argc, char *argv[]) {
    FILE       *fin, *fout;
    int         state, next_state, i, j, k, ind, score, verbose;
    float       symbol;
    uint8_t     bit, bit_buffer[UW_BITS];
    double      symbol_buf[SYMBOLS_PER_PACKET];
    double      symbol_buf_no_rs232[SYMBOLS_PER_PACKET];
    double      llr[SYMBOLS_PER_PACKET];
    char        unpacked_packet[CODELENGTH];
    uint8_t     packet[BYTES_PER_PACKET+CRC_BYTES];
    uint8_t     abyte;
    uint16_t    tx_checksum, rx_checksum, packet_errors, packets;
    int         CodeLength, iter, parityCheckCount;
    struct LDPC ldpc;

    assert(sizeof(uw) == UW_BITS);

    /* LDPC parameters */

    CodeLength = CODELENGTH;                    /* length of entire codeword in bits */

    /* set up LDPC code from include file constants */

    ldpc.max_iter = MAX_ITER;
    ldpc.dec_type = 0;
    ldpc.q_scale_factor = 1;
    ldpc.r_scale_factor = 1;
    ldpc.CodeLength = CODELENGTH;
    ldpc.NumberParityBits = NUMBERPARITYBITS;
    ldpc.NumberRowsHcols = NUMBERROWSHCOLS;
    ldpc.max_row_weight = MAX_ROW_WEIGHT;
    ldpc.max_col_weight = MAX_COL_WEIGHT;
    ldpc.H_rows = H_rows;
    ldpc.H_cols = H_cols;

    /* process command line ----------------------------------------------*/

    if (argc < 3) {
	fprintf(stderr, "usage: drs232 InputOneSymbolPerFloat OutputPackets [-v[v]]\n");
 	exit(1);
    }

    if (strcmp(argv[1], "-")  == 0) fin = stdin;
    else if ( (fin = fopen(argv[1],"rb")) == NULL ) {
	fprintf(stderr, "Error opening input file: %s: %s.\n",
         argv[1], strerror(errno));
	exit(1);
    }

    if (strcmp(argv[2], "-") == 0) fout = stdout;
    else if ( (fout = fopen(argv[2],"wb")) == NULL ) {
	fprintf(stderr, "Error opening output file: %s: %s.\n",
         argv[2], strerror(errno));
	exit(1);
    }

    verbose = 0;
    if (argc > 3) {
        if (strcmp(argv[3], "-v") == 0) {
            verbose = 1;
        }
        if (strcmp(argv[3], "-vv") == 0) {
            verbose = 2;
        }
    }

    state = LOOK_FOR_UW;
    memset(bit_buffer,0,  sizeof(bit_buffer));

    packet_errors = packets = 0;

    while(fread(&symbol, sizeof(float), 1, fin) == 1) {

        /* make hard decision for purpose of UW detection */

        bit = symbol < 0;
        //printf("symbol; %f bit: %d\n", symbol, bit);
        next_state = state;
        if (state == LOOK_FOR_UW) {

            /* put latest input bit into sliding buffer */
            
            for(i=0; i<UW_BITS-1; i++) {
                bit_buffer[i] = bit_buffer[i+1];
            }
            bit_buffer[i] = bit;

            /* lets see if it matches the UW */

            score = 0;
            for(i=0; i<UW_BITS; i++) {
                score += (bit_buffer[i] == uw[i]);
                /* if (i == BITS_PER_BYTE)
                    printf(" ");
                    printf("%1d", unpacked_packet[i]); */
            }
            //printf("\n");
            
            //fprintf(stderr,"UW score: %d\n", score);
            if (score >= (UW_BITS-UW_ALLOWED_ERRORS)) {
                //fprintf(stderr,"UW found! score: %d\n verbose: %d\n", score, verbose);
                ind = 0;
                next_state = COLLECT_PACKET;
            }             
        }

        if (state == COLLECT_PACKET) {
            symbol_buf[ind++] = symbol;
 
            if (ind == SYMBOLS_PER_PACKET) {

               /* OK we have enough bits, remove RS232 sync symbols.
                  This is set up for bit<->byte ordering as per python
                  tx code */

               for(i=0,k=0; i<SYMBOLS_PER_PACKET; i+=BITS_PER_BYTE) {
                   for(j=0; j<8; j++) {
                       symbol_buf_no_rs232[k+j] = symbol_buf[i+7-j+1];
                   }
                   k += 8;
               }

               /* now LDPC decode */

               sd_to_llr(llr, symbol_buf_no_rs232, CodeLength);
               iter = run_ldpc_decoder(&ldpc, unpacked_packet, llr, &parityCheckCount);

               /* pack into bytes */

               for(i=0; i<BYTES_PER_PACKET+CRC_BYTES; i++) {
                   abyte = 0;
                   for(j=0; j<8; j++)
                       abyte |= unpacked_packet[8*i+j] << (7-j);
                   packet[i] = abyte;
               }

               /* then output if CRC check is OK */

               rx_checksum = gen_crc16(packet, BYTES_PER_PACKET);
               tx_checksum = packet[BYTES_PER_PACKET] + (packet[BYTES_PER_PACKET+1] << 8);

               if (verbose == 2) {
                   if (rx_checksum != tx_checksum) {
                       fprintf(stderr, "tx_checksum: 0x%02x rx_checksum: 0x%02x\n", 
                               tx_checksum, rx_checksum);
                   }
               }
               
               packets++;
               if (rx_checksum == tx_checksum) {
                   fwrite(packet, sizeof(char), BYTES_PER_PACKET, fout);
               }
               else
                   packet_errors++;

               if (verbose) {
                   fprintf(stderr, "packets: %d packet_errors: %d PER: %4.3f iter: %d\n", 
                           packets, packet_errors, 
                           (float)packet_errors/packets, iter);
               }
               //exit(0);
               next_state = LOOK_FOR_UW;
            }

        }
        //if (bits_read == (16*10 + UNPACKED_PACKET_BYTES))
        //    exit(0);

        state = next_state;      
    }

    fclose(fin);
    fclose(fout);

    fprintf(stderr, "packets: %d packet_errors: %d PER: %4.3f\n", packets, packet_errors, 
            (float)packet_errors/packets);

    return 0;
}


