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
#include "fsk.h"

#define TEST_FRAME_SIZE 100  /* must match fsk_get_test_bits.c */

int main(int argc,char *argv[]){
    int bitcnt,biterr,i,errs;
    FILE *fin;
    uint8_t *bitbuf_tx, *bitbuf_rx, abit;
    
    if(argc != 2){
        fprintf(stderr,"usage: %s InputBitsOnePerByte \n",argv[0]);
        exit(1);
    }
    
    if(strcmp(argv[1],"-")==0 || argc<2){
        fin = stdin;
    }else{
        fin = fopen(argv[1],"r");
    }
    
    if(fin==NULL){
        fprintf(stderr,"Couldn't open input file: %s\n", argv[1]);
        goto cleanup;
    }

    /* allocate buffers for processing */
    bitbuf_tx = (uint8_t*)alloca(sizeof(uint8_t)*TEST_FRAME_SIZE);
    bitbuf_rx = (uint8_t*)alloca(sizeof(uint8_t)*TEST_FRAME_SIZE);
    
    /* Generate known tx frame from known seed */
    srand(158324);
    for(i=0; i<TEST_FRAME_SIZE; i++){
	bitbuf_tx[i] = rand()&0x1;
	bitbuf_rx[i] = 0;
    }
    
    bitcnt = 0;
    biterr = 0;
    
    while(fread(&abit,sizeof(uint8_t),1,fin)>0){

        /* update silding window of input bits */

        for(i=0; i<TEST_FRAME_SIZE-1; i++) {
            bitbuf_rx[i] = bitbuf_rx[i+1];
        }
        bitbuf_rx[TEST_FRAME_SIZE-1] = abit;

        /* compare to know tx frame.  If they are time aligned, there
           will be a fairly low bit error rate */

        errs = 0;
        for(i=0; i<TEST_FRAME_SIZE; i++) {
            if (bitbuf_rx[i] != bitbuf_tx[i]) {
                errs++;
            }
        }
        if (errs < 0.1*TEST_FRAME_SIZE) {
            /* OK, we have a valid test frame sync, so lets count errors */
            bitcnt += TEST_FRAME_SIZE;
            biterr += errs;
            fprintf(stderr,"errs: %d FSK BER %f, bits tested %d, bit errors %d\n",
                    errs, ((float)biterr/(float)bitcnt),bitcnt,biterr);
        }
    }
 
    
 cleanup:
    fclose(fin);

    return 0;
}
