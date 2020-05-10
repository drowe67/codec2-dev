/*---------------------------------------------------------------------------*\

  FILE........: fsk_get_test_bits.c
  AUTHOR......: Brady O'Brien
  DATE CREATED: January 2016

  Generates a pseudorandom sequence of bits for testing of fsk_mod and
  fsk_demod.

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

#include <stdio.h>
#include <string.h>
#include <getopt.h>
#include "fsk.h"

#define TEST_FRAME_SIZE 100  /* must match fsk_get_test_bits.c */

#define FRAME_DETECTION_THRESHOLD 0.1

int main(int argc,char *argv[]){
    int bitcnt,biterr,i,errs;
    int framesize = TEST_FRAME_SIZE;
    float threshold = FRAME_DETECTION_THRESHOLD;
    float ber_thresh = 1.0;
    FILE *fin;
    uint8_t *bitbuf_tx, *bitbuf_rx, abit;
    int verbose = 1;
    
    char usage[] = "usage: %s InputBitsOnePerByte [-f frameSizeBits] [-t VaildFrameBERThreshold] [-b berPassThreshold] InputFile\n";

    int opt;
    while ((opt = getopt(argc, argv, "f:t:b:hq")) != -1) {
        switch (opt) {
        case 'b':
            ber_thresh = atof(optarg);
            break;
        case 'f':
            framesize = atoi(optarg);
            break;
        case 't':
            threshold = atof(optarg);
            break;
        case 'q':
            verbose = 0;
            break;
        case 'h':
        default:
            fprintf(stderr, usage, argv[0]);
            exit(1);
        }
    }

    char *fname = argv[optind++];
    if ((strcmp(fname,"-")==0) || (argc<2)){
        fin = stdin;
    } else {
        fin = fopen(fname,"r");
    }
    
    if(fin==NULL){
        fprintf(stderr,"Couldn't open input file: %s\n", argv[1]);
        exit(1);
    }

    /* allocate buffers for processing */
    bitbuf_tx = (uint8_t*)alloca(sizeof(uint8_t)*framesize);
    bitbuf_rx = (uint8_t*)alloca(sizeof(uint8_t)*framesize);
    
    /* Generate known tx frame from known seed */
    srand(158324);
    for(i=0; i<framesize; i++){
	bitbuf_tx[i] = rand()&0x1;
	bitbuf_rx[i] = 0;
    }
    
    bitcnt = 0;
    biterr = 0;
    float ber = 0.5;
    
    while(fread(&abit,sizeof(uint8_t),1,fin)>0){

        /* update silding window of input bits */

        for(i=0; i<framesize-1; i++) {
            bitbuf_rx[i] = bitbuf_rx[i+1];
        }
        bitbuf_rx[framesize-1] = abit;

        /* compare to know tx frame.  If they are time aligned, there
           will be a fairly low bit error rate */

        errs = 0;
        for(i=0; i<framesize; i++) {
            if (bitbuf_rx[i] != bitbuf_tx[i]) {
                errs++;
            }
        }
        if (errs < threshold * framesize) {
            /* OK, we have a valid test frame sync, so lets count errors */
            bitcnt += framesize;
            biterr += errs;
            ber =  (float)biterr/(float)bitcnt;
            if (verbose) {
                fprintf(stderr,"errs: %d FSK BER %f, bits tested %d, bit errors %d\n",
                    errs, ber, bitcnt, biterr);
            }
        }
    }
 
    fclose(fin);

    fprintf(stderr,"errs: %d FSK BER %f, bits tested %d, bit errors %d ", errs, ber, bitcnt, biterr);
    if (ber < ber_thresh) {
        fprintf(stderr,"PASS\n");
        return 0;
    }
    else {
        fprintf(stderr,"FAIL\n");
        return 1;
    }
}
