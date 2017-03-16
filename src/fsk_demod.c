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
#include <stdlib.h>
#include <getopt.h>

#include "fsk.h"
#include "codec2_fdmdv.h"
#include "modem_stats.h"

int main(int argc,char *argv[]){
    struct FSK *fsk;
    struct MODEM_STATS stats;
    int Fs,Rs,M,P,stats_ctr,stats_loop;
    float loop_time;
    int enable_stats = 0;
    int hbr = 1;
    FILE *fin,*fout;
    uint8_t *bitbuf;
    int16_t *rawbuf;
    COMP *modbuf;
    float *sdbuf;
    int i,j,Ndft;
    int soft_dec_mode = 0;
    stats_loop = 0;
    int complex_input = 1, bytes_per_sample = 2;
    int stats_rate = 8;
    P = 0;
    M = 0;
    
    /* -m 2/4 --mode 2/4  - FSK mode
     * -l --lbr      - Low bit rate mode
     * -p n --conv n    - Downconverted symbol size
     * If p is unspecified, it will default to Ts
     * -c --cs16      - Complex (signed 16-bit)
     * -d --cu8      - Complex (unsigned 8-bit)
     * If neither -c or -d are specified, input will be real valued signed 16-bit
     * -t --stats - dump demod statistics to stderr
     * -s --soft-dec - ouput soft decision (float 32-bit)
     */
     /* usage: [-l] [-p P]  [-s] [(-c|-d)] [-t] (2|4) SampleRate SymbolRatez InputModemRawFile OutputOneBitPerCharFile */
     
     int o = 0;
     int opt_idx = 0;
     while( o != -1 ){
        static struct option long_opts[] = {
            {"help",      no_argument,        0, 'h'},
            {"lbr",       no_argument,        0, 'l'},
            {"conv",      required_argument,  0, 'p'},
            {"cs16",      no_argument,        0, 'c'},
            {"cu8",       no_argument,        0, 'd'},
            {"stats",     optional_argument,  0, 't'},
            {"soft-dec",  no_argument,        0, 's'},
            {0, 0, 0, 0}
        };
        
        o = getopt_long(argc,argv,"hlp:cdt::s",long_opts,&opt_idx);
        
        switch(o){
            case 'l':
                hbr = 0;
                break;
            case 'c':
                complex_input = 2;
                bytes_per_sample = 2;
                break;
            case 'd':
                complex_input = 2;
                bytes_per_sample = 1;
                break;
            case 't':
                enable_stats = 1;
                if(optarg != NULL){
                    stats_rate = atoi(optarg);
                    if(stats_rate == 0){
                        stats_rate = 8;
                    }
                }
                break;
            case 's':
                soft_dec_mode = 1;
                break;
            case 'p':
                P = atoi(optarg);
                break;
            case 'h':
            case '?':
                goto helpmsg;
                break;
        }
    }
    int dx = optind;
    
    if( (argc - dx) < 5){
        fprintf(stderr, "Too few arguments\n");
        goto helpmsg;
    }
    
    if( (argc - dx) > 5){
        fprintf(stderr, "Too many arguments\n");
        helpmsg:
        fprintf(stderr,"usage: %s [-l] [-p P]  [-s] [(-c|-d)] [-t [r]] (2|4) SampleRate SymbolRate InputModemRawFile OutputFile\n",argv[0]);
        fprintf(stderr," -lP --conv=P      -  P specifies the rate at which symbols are down-converted before further processing\n");
        fprintf(stderr,"                        P must be divisible by the symbol size. Smaller P values will result in faster\n");
        fprintf(stderr,"                        processing but lower demodulation preformance. If no P value is specified,\n");
        fprintf(stderr,"                        P will default to it's highes possible value\n");
        fprintf(stderr," -c --cs16         -  The raw input file will be in complex signed 16 bit format.\n");
        fprintf(stderr," -d --cu8          -  The raw input file will be in complex unsigned 8 bit format.\n");
        fprintf(stderr,"                        If neither -c nor -d are used, the input should be in signed 16 bit format.\n");
        fprintf(stderr," -t[r] --stats=[r] -  Print out modem statistics to stderr in JSON.\n");
        fprintf(stderr,"                         r, if provided, sets the number of modem frames between statistic printouts\n");
        fprintf(stderr," -s --soft-dec     -  The output file will be in a soft-decision format, with one 32-bit float per bit.\n");
        fprintf(stderr,"                        If -s is not used, the output will be in a 1 byte-per-bit format.\n");
        exit(1);
    }
    
    /* Extract parameters */
    M = atoi(argv[dx]);
    Fs = atoi(argv[dx + 1]);
    Rs = atoi(argv[dx + 2]);
    
    if( P == 0 ){
        P = Fs/Rs;
    }
    
    if( (M!=2) && (M!=4) ){
        fprintf(stderr,"Mode %d is not valid. Mode must be 2 or 4.\n",M);
        goto helpmsg;
    }
    
    /* Open files */
    if(strcmp(argv[dx + 3],"-")==0){
        fin = stdin;
    }else{
        fin = fopen(argv[dx + 3],"r");
    }
    
    if(strcmp(argv[dx + 4],"-")==0){
        fout = stdout;
    }else{
        fout = fopen(argv[dx + 4],"w");
    }
    
    
    /* set up FSK */
    if(!hbr)
    fsk = fsk_create(Fs,Rs,M,1200,400);
    else
    fsk = fsk_create_hbr(Fs,Rs,P,M,1200,400);
    
    if(fin==NULL || fout==NULL || fsk==NULL){
        fprintf(stderr,"Couldn't open test vector files\n");
        exit(1);
    }
    
    if(enable_stats){
        loop_time = ((float)fsk_nin(fsk))/((float)Fs);
        stats_loop = (int)(1/(stats_rate*loop_time));
        stats_ctr = 0;
    }
    
    /* allocate buffers for processing */
    if(soft_dec_mode){
        sdbuf = (float*)malloc(sizeof(float)*fsk->Nbits);
    }else{
        bitbuf = (uint8_t*)malloc(sizeof(uint8_t)*fsk->Nbits);
    }
    rawbuf = (int16_t*)malloc(bytes_per_sample*(fsk->N+fsk->Ts*2)*complex_input);
    modbuf = (COMP*)malloc(sizeof(COMP)*(fsk->N+fsk->Ts*2));
    
    /* Demodulate! */
    while( fread(rawbuf,bytes_per_sample*complex_input,fsk_nin(fsk),fin) == fsk_nin(fsk) ){
 
        /* convert input to a buffer of floats.  Note scaling isn't really necessary for FSK */

       if (complex_input == 1) {
            /* S16 real input */
            for(i=0;i<fsk_nin(fsk);i++){
                modbuf[i].real = ((float)rawbuf[i])/FDMDV_SCALE;
                modbuf[i].imag = 0.0;
            }
        }
        else {
            if (bytes_per_sample == 1) {
                /* U8 complex */
                uint8_t *rawbuf_u8 = (uint8_t*)rawbuf;
                for(i=0;i<fsk_nin(fsk);i++){
                    modbuf[i].real = ((float)rawbuf_u8[2*i]-127.0)/128.0;
                    modbuf[i].imag = ((float)rawbuf_u8[2*i+1]-127.0)/128.0;
                }
            }
            else {
                /* S16 complex */
                for(i=0;i<fsk_nin(fsk);i++){
                    modbuf[i].real = ((float)rawbuf[2*i])/FDMDV_SCALE;
                    modbuf[i].imag = ((float)rawbuf[2*i+1]/FDMDV_SCALE);
                }
            }            
        }
        if(soft_dec_mode){
            fsk_demod_sd(fsk,sdbuf,modbuf);
        }else{
            fsk_demod(fsk,bitbuf,modbuf);
        }
        
        if(enable_stats && stats_ctr <= 0){
            fsk_get_demod_stats(fsk,&stats);
	    /* Print standard 2FSK stats */
            fprintf(stderr,"{\"EbNodB\": %2.2f,\t\"ppm\": %d,",stats.snr_est,(int)fsk->ppm);
            fprintf(stderr,"\t\"f1_est\":%.1f,\t\"f2_est\":%.1f",fsk->f_est[0],fsk->f_est[1]);
	    /* Print 4FSK stats if in 4FSK mode */
            if(fsk->mode == 4){
                fprintf(stderr,",\t\"f3_est\":%.1f,\t\"f4_est\":%.1f",fsk->f_est[2],fsk->f_est[3]);
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
    
    fclose(fin);
    fclose(fout);
    fsk_destroy(fsk);

    return 0;
}

