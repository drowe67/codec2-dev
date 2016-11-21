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

#include "codec2_fdmdv.h"
#include "modem_stats.h"

int main(int argc,char *argv[]){
    struct FSK *fsk;
    struct MODEM_STATS stats;
    int Fs,Rs,M,P,stats_ctr,stats_loop;
    float loop_time;
    int enable_stats = 0;
    int hbr = 0;
    FILE *fin,*fout;
    uint8_t *bitbuf;
    int16_t *rawbuf;
    COMP *modbuf;
    float *sdbuf;
    int i,j,Ndft;
    int soft_dec_mode = 0;
    stats_loop = 0;
    int complex_input = 1;
    
    if(argc<7){
        fprintf(stderr,"usage: %s (2|4|2X|4X|2XS) P SampleFreq SymbolFreq InputModemRawFile OutputOneBitPerCharFile [S] [C]\n",argv[0]);
        fprintf(stderr, "2   - 2FSK low bit rate\n4   - 4FSK low bit rate\n2X  - 2FSK high bit rate\n4   - 4FSK high bit rate\n2XS - high bit rate soft decision output\n\n");
         fprintf(stderr, "P   - timing estimator window size, see README_fsk\n");
         fprintf(stderr, "S   - dump demod stats to stderr for plotting with octave/fskdemodgui.py\n");
         fprintf(stderr, "C   - complex (two sample) input\n");
       exit(1);
    }
    
    /* Extract parameters */
    P  = atoi(argv[2]);
    Fs = atoi(argv[3]);
    Rs = atoi(argv[4]);
    
    /* Open files */
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

    /* Handle high-bit-rate special cases */
    if(strcmp(argv[1],"2X")==0){
        M = 2;
        hbr = 1;
    }else if(strcmp(argv[1],"2XS")==0){
        M = 2;
        hbr = 1;
        soft_dec_mode = 1;
    }else if(strcmp(argv[1],"4X")==0){
        M = 4;
        hbr = 1;
    }else {
        M = atoi(argv[1]);
    }
    
    /* set up FSK */
    if(!hbr)
    fsk = fsk_create(Fs,Rs,M,1200,400);
    else
    fsk = fsk_create_hbr(Fs,Rs,P,M,1200,400);
    
    if(fin==NULL || fout==NULL || fsk==NULL){
        fprintf(stderr,"Couldn't open test vector files\n");
        goto cleanup;
    }
    
    /* Check for and enable stat printing and complex input */
    /* TODO: design better command line arguments */
    if(argc>7){
        if(strcmp(argv[7],"S")==0){
            enable_stats = 1;
            fsk_setup_modem_stats(fsk,&stats);
            loop_time = ((float)fsk_nin(fsk))/((float)Fs);
            stats_loop = (int)(.125/loop_time);
            stats_ctr = 0;
        }
        if(strcmp(argv[7],"C")==0){
            complex_input = 2;
        }
    }
    if(argc>8){
        if(strcmp(argv[8],"S")==0){
            enable_stats = 1;
            fsk_setup_modem_stats(fsk,&stats);
            loop_time = ((float)fsk_nin(fsk))/((float)Fs);
            stats_loop = (int)(.125/loop_time);
            stats_ctr = 0;
        }
        if(strcmp(argv[8],"C")==0){
            complex_input = 2;
        }
    }
       
    /* allocate buffers for processing */
    if(soft_dec_mode){
        sdbuf = (float*)malloc(sizeof(float)*fsk->Nbits);
    }else{
        bitbuf = (uint8_t*)malloc(sizeof(uint8_t)*fsk->Nbits);
    }
    rawbuf = (int16_t*)malloc(sizeof(int16_t)*(fsk->N+fsk->Ts*2)*complex_input);
    modbuf = (COMP*)malloc(sizeof(COMP)*(fsk->N+fsk->Ts*2));
    
    /* Demodulate! */
    while( fread(rawbuf,sizeof(int16_t)*complex_input,fsk_nin(fsk),fin) == fsk_nin(fsk) ){
        if (complex_input == 1) {
            for(i=0;i<fsk_nin(fsk);i++){
                modbuf[i].real = ((float)rawbuf[i])/FDMDV_SCALE;
                modbuf[i].imag = 0.0;
            }
        }
        else {
            for(i=0;i<fsk_nin(fsk);i++){
                modbuf[i].real = ((float)rawbuf[2*i])/FDMDV_SCALE;
                modbuf[i].imag = ((float)rawbuf[2*i+1])/FDMDV_SCALE;
            }
        }
        if(soft_dec_mode){
            fsk_demod_sd(fsk,sdbuf,modbuf);
        }else{
            fsk_demod(fsk,bitbuf,modbuf);
        }
        
        if(enable_stats && stats_ctr <= 0){
	    /* Print standard 2FSK stats */
            fprintf(stderr,"{\"EbNodB\": %2.2f,\t\"ppm\": %d,",stats.snr_est,(int)fsk->ppm);
            fprintf(stderr,"\t\"f1_est\":%.1f,\t\"f2_est\":%.1f",fsk->f1_est,fsk->f2_est);
	    /* Print 4FSK stats if in 4FSK mode */
            if(fsk->mode == 4){
                fprintf(stderr,",\t\"f3_est\":%.1f,\t\"f4_est\":%.1f",fsk->f3_est,fsk->f4_est);
            }
	    
	    /* Print the eye diagram */
            fprintf(stderr,",\t\"eye_diagram\":[");
            for(i=0;i<stats.neyetr;i++){
                fprintf(stderr,"[");
                for(j=0;j<stats.neyesamp;j++){
                    fprintf(stderr,"%f ",stats.rx_eye[i][j]);
                    if(j<stats.neyesamp-1) fprintf(stderr,",");
                }
                fprintf(stderr,"]");
                if(i<stats.neyetr-1) fprintf(stderr,",");
            }
            fprintf(stderr,"],");
	    
	    /* Print a sample of the FFT from the freq estimator */
	    fprintf(stderr,"\"samp_fft\":[");
	    Ndft = fsk->Ndft/2;
	    for(i=0; i<Ndft; i++){
		fprintf(stderr,"%f ",(fsk->fft_est)[i]);
		if(i<Ndft-1) fprintf(stderr,",");
	    }
	    fprintf(stderr,"]}\n");
            stats_ctr = stats_loop;
        }
        stats_ctr--;
            /*for(i=0;i<fsk->Nbits;i++){
            t = (int)bitbuf[i];
        }*/
        if(soft_dec_mode){
            fwrite(sdbuf,sizeof(float),fsk->Nbits,fout);
        }else{
            fwrite(bitbuf,sizeof(uint8_t),fsk->Nbits,fout);
        }
            if(fin == stdin || fout == stdin){
            fflush(fin);
            fflush(fout);
        }
    }
    
    if(soft_dec_mode){
        free(sdbuf);
    }else{
        free(bitbuf);
    }
    
    free(rawbuf);
    free(modbuf);
    
    cleanup:
    fclose(fin);
    fclose(fout);
    fsk_destroy(fsk);
    exit(0);
}

