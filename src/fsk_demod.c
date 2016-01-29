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
#include "codec2_fdmdv.h"

int main(int argc,char *argv[]){
    struct FSK *fsk;
    int Fs,Rs,M;
    FILE *fin,*fout;
    uint8_t *bitbuf;
    int16_t *rawbuf;
    float *modbuf;
    int i,t;
    
    if(argc<6){
        fprintf(stderr,"usage: %s Mode SampleFreq SymbolFreq InputModemRawFile OutputOneBitPerCharFile [OctaveLogFile]\n",argv[0]);
        exit(1);
    }
    
    /* Extract parameters */
    M  = atoi(argv[1]);
    Fs = atoi(argv[2]);
    Rs = atoi(argv[3]);
    
    /* Open files */
    if(strcmp(argv[4],"-")==0){
		fin = stdin;
	}else{
		fin = fopen(argv[4],"r");
	}
	
	if(strcmp(argv[5],"-")==0){
		fout = stdout;
	}else{
		fout = fopen(argv[5],"w");
	}

    
    if(argc>6)
		modem_probe_init("fsk2",argv[6]);
	
    /* set up FSK */
    fsk = fsk_create(Fs,Rs,M,1200,400);
    
    if(fin==NULL || fout==NULL || fsk==NULL){
        fprintf(stderr,"Couldn't open test vector files\n");
        goto cleanup;
    }
    
    /* allocate buffers for processing */
    bitbuf = (uint8_t*)alloca(sizeof(uint8_t)*fsk->Nbits);
    rawbuf = (int16_t*)alloca(sizeof(int16_t)*(fsk->N+fsk->Ts*2));
    modbuf = (float*)alloca(sizeof(float)*(fsk->N+fsk->Ts*2));
    
    /* Demodulate! */
    while( fread(rawbuf,sizeof(int16_t),fsk_nin(fsk),fin) == fsk_nin(fsk) ){
		for(i=0;i<fsk_nin(fsk);i++){
			modbuf[i] = ((float)rawbuf[i])/FDMDV_SCALE;
		}
		modem_probe_samp_f("t_d_sampin",modbuf,fsk_nin(fsk));
        fsk_demod(fsk,bitbuf,modbuf);
        for(i=0;i<fsk->Nbits;i++){
			t = (int)bitbuf[i];
			modem_probe_samp_i("t_d_bitout",&t,1);
		}
        fwrite(bitbuf,sizeof(uint8_t),fsk->Nbits,fout);
        
        if(fin == stdin || fout == stdin){
			fflush(fin);
			fflush(fout);
		}
    }
    
    modem_probe_close();
    cleanup:
    fclose(fin);
    fclose(fout);
    fsk_destroy(fsk);
    exit(0);
}

