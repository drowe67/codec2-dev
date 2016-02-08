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
#include "fmfsk.h"

#define MODEMPROBE_ENABLE
#include "modem_probe.h"
#include "codec2_fdmdv.h"

int main(int argc,char *argv[]){
    struct FMFSK *fmfsk;
    int Fs,Rb;
    FILE *fin,*fout;
    uint8_t *bitbuf;
    int16_t *rawbuf;
    float *modbuf;
    int i,t;
    
    if(argc<4){
        fprintf(stderr,"usage: %s SampleFreq BitRate InputModemRawFile OutputOneBitPerCharFile [OctaveLogFile]\n",argv[0]);
        exit(1);
    }
    
    /* Extract parameters */
    Fs = atoi(argv[1]);
    Rb = atoi(argv[2]);
    
    /* Open files */
    if(strcmp(argv[3],"-")==0){
	fin = stdin;
    }else{
	fin = fopen(argv[3],"r");
    }
	
    if(strcmp(argv[4],"-")==0){
	fout = stdout;
    }else{
	fout = fopen(argv[4],"w");
    }

    
    if(argc>4)
	modem_probe_init("fmfsk2",argv[5]);
	
    /* set up FSK */
    fmfsk = fmfsk_create(Fs,Rb);
    
    if(fin==NULL || fout==NULL || fmfsk==NULL){
        fprintf(stderr,"Couldn't open test vector files\n");
        goto cleanup;
    }
    
    /* allocate buffers for processing */
    bitbuf = (uint8_t*)alloca(sizeof(uint8_t)*fmfsk->nbit);
    rawbuf = (int16_t*)alloca(sizeof(int16_t)*(fmfsk->N+fmfsk->Ts*2));
    modbuf = (float*)alloca(sizeof(float)*(fmfsk->N+fmfsk->Ts*2));
    
    /* Demodulate! */
    while( fread(rawbuf,sizeof(int16_t),fmfsk_nin(fmfsk),fin) == fmfsk_nin(fmfsk) ){
	for(i=0;i<fmfsk_nin(fmfsk);i++){
	    modbuf[i] = ((float)rawbuf[i])/FDMDV_SCALE;
	}
	
	modem_probe_samp_f("t_d_sampin",modbuf,fmfsk_nin(fmfsk));
        fmfsk_demod(fmfsk,bitbuf,modbuf);
        
	for(i=0;i<fmfsk->nbit;i++){
	    t = (int)bitbuf[i];
	    modem_probe_samp_i("t_d_bitout",&t,1);
	}
        
	fwrite(bitbuf,sizeof(uint8_t),fmfsk->nbit,fout);
        
        if(fin == stdin || fout == stdin){
	    fflush(fin);
	    fflush(fout);
	}
    }
    
    modem_probe_close();
    cleanup:
    fclose(fin);
    fclose(fout);
    fmfsk_destroy(fmfsk);
    exit(0);
}

