/*---------------------------------------------------------------------------*\

  FILE........: fsk_demod.c
  AUTHOR......: Brady O'Brien
  DATE CREATED: 8 January 2016

  C test driver for fsk_demod in fsk.c. Reads in a stream of 32 bit cpu endian
  floats and writes out the detected bits
   

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
#include "fsk.h"

#define MODEMPROBE_ENABLE
#include "modem_probe.h"

int main(int argc,char *argv[]){
    struct FSK *fsk;
    int Fs,Rs,f1,f2;
    FILE *fin,*fout;
    uint8_t *bitbuf;
    float *modbuf;
    
    if(argc<5){
        printf("usage: %s SampleFreq SymbolFreq InputModemRawFile OutputOneBitPerCharFile [OctaveLogFile]\n",argv[0]);
        exit(1);
    }
    
    /* Extract parameters */
    Fs = atoi(argv[1]);
    Rs = atoi(argv[2]);
    
    /* Open files */
    fin = fopen(argv[3],"r");
    fout = fopen(argv[4],"w");
    
    if(argc>5)
		modem_probe_init("fsk2",argv[5]);
	
    /* set up FSK */
    fsk = fsk_create(Fs,Rs,0,0);
    
    if(fin==NULL || fout==NULL || fsk==NULL){
        printf("Couldn't open test vector files\n");
        goto cleanup;
    }
    
    /* allocate buffers for processing */
    bitbuf = (uint8_t*)alloca(sizeof(uint8_t)*fsk->Nsym);
    modbuf = (float*)alloca(sizeof(float)*(fsk->N+fsk->Ts*2));
    
    /* Demodulate! */
    while( fread(modbuf,sizeof(float),fsk_nin(fsk),fin) == fsk_nin(fsk) ){
        fsk_demod(fsk,bitbuf,modbuf);
        fwrite(bitbuf,sizeof(uint8_t),fsk->Nsym,fout);
    }
    
    modem_probe_close();
    cleanup:
    fclose(fin);
    fclose(fout);
    fsk_destroy(fsk);
    exit(0);
}

