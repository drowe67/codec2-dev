/*---------------------------------------------------------------------------*\

  FILE........: fsk_mod.c
  AUTHOR......: Brady O'Brien
  DATE CREATED: 8 January 2016

  C test driver for fsk_mod in fsk.c. Reads in a set of bits to modulate
   from a file, passed as a parameter, and writes modulated output to
   another file
   

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
#include "codec2_fdmdv.h"

int main(int argc,char *argv[]){
    struct FSK *fsk;
    int Fs,Rs,f1,f2;
    int i;
    FILE *fin,*fout;
    uint8_t *bitbuf;
    int16_t *rawbuf;
    float *modbuf;
    
    if(argc<7){
        fprintf(stderr,"usage: %s SampleFreq SymbolFreq TxFreq1 TxFreq2 InputOneBitPerCharFile OutputModRawFile\n",argv[0]);
        exit(1);
    }
    
    /* Extract parameters */
    Fs = atoi(argv[1]);
    Rs = atoi(argv[2]);
    f1 = atoi(argv[3]);
    f2 = atoi(argv[4]);
    
    if(strcmp(argv[5],"-")==0){
		fin = stdin;
	}else{
		fin = fopen(argv[5],"r");
	}
	
	if(strcmp(argv[6],"-")==0){
		fout = stdout;
	}else{
		fout = fopen(argv[6],"w");
	}
    
    
    /* set up FSK */
    fsk = fsk_create(Fs,Rs,f1,f2);
    
    if(fin==NULL || fout==NULL || fsk==NULL){
        fprintf(stderr,"Couldn't open test vector files\n");
        goto cleanup;
    }
    
    /* allocate buffers for processing */
    bitbuf = (uint8_t*)alloca(sizeof(uint8_t)*fsk->Nsym);
    rawbuf = (int16_t*)alloca(sizeof(int16_t)*fsk->N);
    modbuf = (float*)alloca(sizeof(float)*fsk->N);
    
    /* Modulate! */
    while( fread(bitbuf,sizeof(uint8_t),fsk->Nsym,fin) == fsk->Nsym ){
        fsk_mod(fsk,modbuf,bitbuf);
        for(i=0; i<fsk->N; i++){
			rawbuf[i] = (int16_t)(modbuf[i]*(float)FDMDV_SCALE);
		}
        fwrite(rawbuf,sizeof(int16_t),fsk->N,fout);
        
		if(fin == stdin || fout == stdin){
			fflush(fin);
			fflush(fout);
		}
    }
    
    cleanup:
    fclose(fin);
    fclose(fout);
    fsk_destroy(fsk);
    exit(0);
}
