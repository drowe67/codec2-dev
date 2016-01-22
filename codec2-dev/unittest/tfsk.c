/*---------------------------------------------------------------------------*\

  FILE........: tfsk.c
  AUTHOR......: Brady O'Brien
  DATE CREATED: 20 January 2016

  C test driver for fsk_mod and fsk_demod in fsk.c. Reads a file with input
  bits/rf and spits out modulated/demoduladed samples and a dump of internal
  state. To run unit test, see octave/tfsk.m

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


#define MODEMPROBE_ENABLE

#include "modem_probe.h"
#include <stdio.h>

/* Note: This is a dirty hack to force fsk.c to compile with modem probing enabled */
#include "fsk.c"


int main(int argc,char *argv[]){
    struct FSK *fsk;
    int Fs,Rs,f1,f2;
    FILE *fin,*fout;
    uint8_t *bitbuf;
    float *modbuf;
    
    if(argc<7){
        printf("Usage: %s samplerate bitrate f1 f2 infile outfile [probefile]\n",argv[0]);
        exit(1);
    }
    
    /* Extract parameters */
    Fs = atoi(argv[1]);
    Rs = atoi(argv[2]);
    f1 = atoi(argv[3]);
    f2 = atoi(argv[4]);
    
    /* Open files */
    fin = fopen(argv[5],"r");
    fout = fopen(argv[6],"w");
    
    if(argc>7)
		modem_probe_init("fsk2",argv[7]);
	
    /* set up FSK */
    fsk = fsk_create(Fs,Rs,f1,f2);
    
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

