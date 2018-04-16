/*---------------------------------------------------------------------------*\

  FILE........: horus_demod.c
  AUTHOR......: David Rowe
  DATE CREATED: April 2018

  Command line demo program for the Horus API, that exercises
  horus_api.c using file input/output (can be stdin/stdout for real
  time operation).  Prints JSON stats, just like Brady's fsk_demod.c
 
  Can operate in Horus RTTY or Binary mode.

  Testing with a 8000Hz sample rate wave file:

    $ sox ~/Desktop/horus.wav -r 48000 -t raw - | ./horus_demod -m RTTY -v - /dev/nul/

    $ sox ~/Desktop/4FSK_binary_100Rb_8khzfs.wav -r 48000 -t raw - | ./horus_demod -m binary  - -

\*---------------------------------------------------------------------------*/

/*
  Copyright (C) 2018 David Rowe

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
#include <stdio.h>
#include <getopt.h>

#include "horus_api.h"
#include "fsk.h"
#include "horus_l2.h"

int main(int argc, char *argv[]) {
    struct   horus *hstates;
    struct   MODEM_STATS stats;
    FILE    *fin,*fout;
    int      i,j,Ndft,mode;
    int      stats_ctr,stats_loop, stats_rate, verbose, crc_results;
    float    loop_time;
    int      enable_stats = 0;

    stats_loop = 0;
    stats_rate = 8;
    mode = -1;
    verbose = crc_results = 0;
    
    int o = 0;
    int opt_idx = 0;
    while ( o != -1 ) {
        static struct option long_opts[] = {
            {"help",      no_argument,        0, 'h'},
            {"mode",      required_argument,  0, 'm'},
            {"stats",     optional_argument,  0, 't'},
            {0, 0, 0, 0}
        };
        
        o = getopt_long(argc,argv,"hvcm:t::",long_opts,&opt_idx);
        
        switch(o) {
            case 'm':
                if ((strcmp(optarg, "RTTY") == 0) || (strcmp(optarg, "rtty") == 0)) {
                    mode = HORUS_MODE_RTTY;
                }
                if ((strcmp(optarg, "BINARY") == 0) || (strcmp(optarg, "binary") == 0)) {
                    mode = HORUS_MODE_BINARY;
                }
                if (mode == -1) {
                    fprintf(stderr, "use --mode RTTY or --mode binary\n");
                    exit(1);
                }
                break;
            case 't':
                enable_stats = 1;
                if (optarg != NULL){
                    stats_rate = atoi(optarg);
                    if (stats_rate == 0) {
                        stats_rate = 8;
                    }
                }
                break;
            case 'v':
                verbose = 1;
            break;    
            case 'c':
                crc_results = 1;
            break;    
            case 'h':
            case '?':
                goto helpmsg;
                break;
        }
    }
    
    int dx = optind;
    
    if( (argc - dx) < 1) {
        fprintf(stderr, "Too few arguments\n");
        goto helpmsg;
    }
    
    if( (argc - dx) > 5){
        fprintf(stderr, "Too many arguments\n");
    helpmsg:
        fprintf(stderr,"usage: %s -m RTTY|binary [-v] [-c] [-t [r]] InputModemRawFile OutputAsciiFile\n",argv[0]);
        fprintf(stderr,"\n");
        fprintf(stderr,"InputModemRawFile      48 kHz 16 bit shorts real modem signal from radio\n");
        fprintf(stderr," -m RTTY|binary\n"); 
        fprintf(stderr,"--mode=RTTY|binary     RTTY or binary Horus protcols\n");
        fprintf(stderr," -t[r] --stats=[r]     Print out modem statistics to stderr in JSON.\n");
        fprintf(stderr,"                       r, if provided, sets the number of modem frames\n"
                       "                       between statistic printouts\n");
        fprintf(stderr," -v                    verbose debug info\n");
        fprintf(stderr," -c                    display CRC results for each packet\n");
        exit(1);
    }
        
    /* Open files */

    if (verbose) {
         fprintf(stderr, "mode: %d verbose: %d stats_loop: %d stats_rate: %d\n",mode, verbose, stats_loop, stats_rate);
    }
    if (strcmp(argv[dx],"-")==0) {
        fin = stdin;
    } else {
        fin = fopen(argv[dx],"rb");
    }
    
    if (strcmp(argv[dx + 1],"-")==0) {
        fout = stdout;
    } else {
        fout = fopen(argv[dx + 1],"w");
    }
        
    if ((fin==NULL) || (fout==NULL)) {
        fprintf(stderr,"Couldn't open test vector files\n");
        exit(1);
    }

    /* end command line processing */

    hstates = horus_open(mode);
    horus_set_verbose(hstates, verbose);
    
    if (hstates == NULL) {
        fprintf(stderr, "Couldn't open Horus API\n");
        exit(1);
    }
    
    if (enable_stats) {
        loop_time = (float)horus_nin(hstates)/horus_get_Fs(hstates);
        stats_loop = (int)(1.0/(stats_rate*loop_time));
        stats_ctr = 0;
    }
    
    int   max_demod_in = horus_get_max_demod_in(hstates);
    short demod_in[max_demod_in];
    int   max_ascii_out = horus_get_max_ascii_out_len(hstates);
    char  ascii_out[max_ascii_out];
    
    /* Main loop ----------------------------------------------------------------------- */

    while(fread(demod_in, sizeof(short), horus_nin(hstates), fin) ==  horus_nin(hstates)) {

        if (verbose) {
            fprintf(stderr, "read nin %d\n", horus_nin(hstates));
        }
        if (horus_rx(hstates, ascii_out, demod_in)) {
            fprintf(stdout, "%s", ascii_out);
            if (crc_results) {
                if (horus_crc_ok(hstates)) {
                    fprintf(stdout, "  CRC OK");
                } else {
                    fprintf(stdout, "  CRC BAD");
                }
            }
            fprintf(stdout, "\n");
        }
        
        if (enable_stats && stats_ctr <= 0) {

            horus_get_modem_extended_stats(hstates, &stats);

	    /* Print standard 2FSK stats */

            fprintf(stderr,"{\"EbNodB\": %2.2f,\t\"ppm\": %d,",stats.snr_est, (int)stats.clock_offset);
            fprintf(stderr,"\t\"f1_est\":%.1f,\t\"f2_est\":%.1f",stats.f_est[0], stats.f_est[1]);

	    /* Print 4FSK stats if in 4FSK mode */

            if (horus_get_mFSK(hstates) == 4) {
                fprintf(stderr,",\t\"f3_est\":%.1f,\t\"f4_est\":%.1f", stats.f_est[2], stats.f_est[3]);
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
	    
	    fprintf(stderr,"\"samp_fft\":[");

            #ifdef FIXME_LATER
            /* TODO: need a horus_ function to dig into modem spectrum */
            
	    /* Print a sample of the FFT from the freq estimator */

	    Ndft = hstates->fsk->Ndft/2;
	    for(i=0; i<Ndft; i++){
		fprintf(stderr,"%f ",(hstates->fsk->fft_est)[i]);
		if(i<Ndft-1) fprintf(stderr,",");
	    }
            #else

            /* All zero dummy data for now */

 	    Ndft = 128;
	    for(i=0; i<Ndft; i++) {
		fprintf(stderr,"%f ", 0.0);
		if(i<Ndft-1) fprintf(stderr,",");
	    }
            
            #endif

	    fprintf(stderr,"]}\n");
            stats_ctr = stats_loop;
        }
        stats_ctr--;

        if (fin == stdin || fout == stdin){
            fflush(fin);
            fflush(fout);
        }
    }
    
    horus_close(hstates);

    return 0;
}

