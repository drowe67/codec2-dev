/*---------------------------------------------------------------------------*\

  FILE........: ofdm_mod.c
  AUTHOR......: David Rowe
  DATE CREATED: March 2018

  Given an input file of bits (note one bit per char format), outputs
  a raw file (8kHz, 16 bit shorts) of OFDM modem samples ready to send
  over a HF radio channel.

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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <errno.h>

#include "codec2_ofdm.h"
#include "ofdm_internal.h"
#include "interldpc.h"
#include "gp_interleaver.h"

extern int payload_data_bits[];
extern int test_bits_ofdm[];

static struct OFDM_CONFIG *ofdm_config;

static int ofdm_bitsperframe;
static int ofdm_nuwbits;
static int ofdm_ntxtbits;

int opt_exists(char *argv[], int argc, char opt[]) {
    int i;
    for (i=0; i<argc; i++) {
        if (strcmp(argv[i], opt) == 0) {
            return i;
        }
    }
    return 0;
}

int main(int argc, char *argv[])
{
    FILE          *fin, *fout;
    struct OFDM   *ofdm;
    struct LDPC   ldpc;
    int           frames;
    int           i, j, arg;
    int           testframes, ldpc_en, interleave_frames;
   
    if (argc < 3) {
        fprintf(stderr, "\n");
	fprintf(stderr, "usage: %s InputOneCharPerBitFile OutputModemRawFile [--ldpc] [--interleaver depth]\n\n", argv[0]);
        fprintf(stderr, "  --nc         Nsecs  Number of Carriers (17 default)\n");
        fprintf(stderr, "  --tcp        Nsecs  Cyclic Prefix Duration (.002 default)\n");
        fprintf(stderr, "  --ts         Nsecs  Symbol Duration (.018 default)\n");
        fprintf(stderr, "  --testframes Nsecs  Transmit test frames (adjusts test frames for raw and LDPC modes)\n");
        fprintf(stderr, "  --ldpc              Run LDPC decoder.  This forces 112, one char/bit output values\n"
                        "                      per frame.  In testframe mode (-t) raw and coded errors will be counted\n");
        fprintf(stderr, "  --interleave depth  Interleave depth for LDPC frames, e.g. 1,2,4,8,16, default is 1\n");
        fprintf(stderr, "  --txbpf             Transmit band pass filter\n");
        fprintf(stderr, "\n");
	exit(1);
    }

    if (strcmp(argv[1], "-")  == 0)
        fin = stdin;
    else if ( (fin = fopen(argv[1],"rb")) == NULL ) {
        fprintf(stderr, "Error opening input file: %s: %s.\n", argv[1], strerror(errno));
        exit(1);
    }

    if (strcmp(argv[2], "-") == 0)
        fout = stdout;
    else if ( (fout = fopen(argv[2],"wb")) == NULL ) {
        fprintf(stderr, "Error opening output modem sample file: %s: %s.\n", argv[2], strerror(errno));
        exit(1);
    }
    
    if ((ofdm_config = (struct OFDM_CONFIG *) calloc(1, sizeof (struct OFDM_CONFIG))) == NULL) {
        fprintf(stderr, "Out of Memory\n");
        exit(1);
    }

    ofdm_config->centre = 1500.0;
    ofdm_config->fs = 8000.0;			/* Sample Frequency */
    ofdm_config->ofdm_timing_mx_thresh = 0.30;
    ofdm_config->ftwindowwidth = 11;
    ofdm_config->state_str = 16; 		/* state string length */
    ofdm_config->bps = 2;   			/* Bits per Symbol */
    ofdm_config->txtbits = 4; 			/* number of auxiliary data bits */
    ofdm_config->ns = 8;  			/* Number of Symbol frames */

    if ((arg = opt_exists(argv, argc, "--nc"))) {
        ofdm_config->nc = atoi(argv[arg+1]);	/* Number of carriers */
    } else {
        ofdm_config->nc = 17;
    }

    if ((arg = opt_exists(argv, argc, "--tcp"))) {
        ofdm_config->tcp = atoi(argv[arg+1]);	/* Cyclic Prefix duration */
    } else {
        ofdm_config->tcp = 0.0020;
    }

    if ((arg = opt_exists(argv, argc, "--ts"))) {
        ofdm_config->ts = atoi(argv[arg+1]);	/* Symbol duration */
    } else {
        ofdm_config->ts = 0.0180;
    }

    ofdm_config->rs = (1.0f / ofdm_config->ts); 	/* Symbol Rate */

    ofdm = ofdm_create(ofdm_config);
    assert(ofdm != NULL);

    free(ofdm_config);

    /* Get a copy of the actual modem config */
    ofdm_config = ofdm_get_config_param();

    ofdm_bitsperframe = ofdm_get_bits_per_frame();
    ofdm_nuwbits = (ofdm_config->ns - 1) * ofdm_config->bps - ofdm_config->txtbits;
    ofdm_ntxtbits = ofdm_config->txtbits;

    /* Set up default LPDC code.  We could add other codes here if we like */
    
    set_up_hra_112_112(&ldpc, ofdm_config);

    int data_bits_per_frame = ldpc.data_bits_per_frame;
    int coded_bits_per_frame = ldpc.coded_bits_per_frame;
    int Nbitsperframe;

    /* set for LDPC coded or uncoded frames */
    
    ldpc_en = 0; interleave_frames = 1;

    if (opt_exists(argv, argc, "--ldpc")) {

        assert((ofdm_nuwbits+ofdm_ntxtbits+coded_bits_per_frame) == ofdm_bitsperframe); /* sanity check */

        ldpc_en = 1;

        if ((arg = opt_exists(argv, argc, "--interleave"))) {
            interleave_frames = atoi(argv[arg+1]);
        }

        Nbitsperframe = interleave_frames*data_bits_per_frame;
    } else {
        /* vanilla uncoded input bits mode */
        Nbitsperframe = ofdm_bitsperframe;
    }
    
    int Nsamperframe = ofdm_get_samples_per_frame();
    fprintf(stderr, "Nbitsperframe: %d interleave_frames: %d\n", Nbitsperframe, interleave_frames);

    uint8_t tx_bits_char[Nbitsperframe];
    short   tx_scaled[Nsamperframe];
    uint8_t txt_bits_char[ofdm_ntxtbits*interleave_frames];

    for(i=0; i< ofdm_ntxtbits*interleave_frames; i++) {
        txt_bits_char[i] = 0;
    }
   
    testframes = 0;
    int Nframes = 0;

    if ((arg = (opt_exists(argv, argc, "--testframes")))) {
        testframes = 1;
        int Nsec, Nrows;
        Nsec = atoi(argv[arg+1]);
        Nrows = Nsec*ofdm_config->rs;
        Nframes = floor((Nrows-1)/ofdm_config->ns);

        fprintf(stderr, "Nframes: %d\n", Nframes);
    }

    if (opt_exists(argv, argc, "--txbpf")) {
        ofdm_set_tx_bpf(ofdm, 1);
    }

    /* main loop ----------------------------------------------------------------*/
    
    frames = 0;

    while(fread(tx_bits_char, sizeof(char), Nbitsperframe, fin) == Nbitsperframe) {
        if (ldpc_en) {
            /* fancy interleaved LDPC encoded frames ----------------------------------------*/
            
            /* optionally overwrite input data with test frame of
               payload data bits known to demodulator */
            
            if (testframes) {
                uint16_t r[data_bits_per_frame];

                ofdm_rand(r, data_bits_per_frame);

                for (j=0; j<interleave_frames; j++) {
                    for(i=0; i<data_bits_per_frame; i++) {
                        //tx_bits_char[j*data_bits_per_frame + i] = payload_data_bits[i];
                        tx_bits_char[j*data_bits_per_frame + i] = r[i]>16384;
                    }
                }
            }

            complex float tx_sams[interleave_frames*Nsamperframe];
            ofdm_ldpc_interleave_tx(ofdm, &ldpc, tx_sams, tx_bits_char, txt_bits_char, interleave_frames, ofdm_config);

            for (j=0; j<interleave_frames; j++) { 
                for(i=0; i<Nsamperframe; i++) {
                    tx_scaled[i] = OFDM_AMP_SCALE * crealf(tx_sams[j*Nsamperframe+i]);
                }

                fwrite(tx_scaled, sizeof(short), Nsamperframe, fout);
                frames++;
            }
         } else {
            /* just modulate uncoded raw bits ----------------------------------------------*/
            
            if (testframes) {
                /* build up a test frame consisting of unique word, txt bits, and psuedo-random
                   uncoded payload bits.  The psuedo-random generator is the same as Octave so
                   it can interoperate with ofdm_tx.m/ofdm_rx.m */

                int Npayloadbits = Nbitsperframe-(ofdm_nuwbits+ofdm_ntxtbits);
                uint16_t r[Npayloadbits];
                uint8_t  payload_bits[Npayloadbits];

                ofdm_rand(r, Npayloadbits);

                for(i=0; i<Npayloadbits; i++) {
                    payload_bits[i] = r[i] > 16384;
                    //fprintf(stderr,"%d %d ", r[j], tx_bits_char[i]);
                }

                uint8_t txt_bits[ofdm_ntxtbits];

                for(i=0; i<ofdm_ntxtbits; i++) {
                    txt_bits[i] = 0;
                }

                ofdm_assemble_modem_frame(tx_bits_char, payload_bits, txt_bits);
            }

            int tx_bits[Nbitsperframe];

            for(i=0; i<Nbitsperframe; i++)
                tx_bits[i] = tx_bits_char[i];

            COMP tx_sams[Nsamperframe];
            ofdm_mod(ofdm, tx_sams, tx_bits);

            /* scale and save to disk as shorts */

            for(i=0; i<Nsamperframe; i++)
                tx_scaled[i] = OFDM_AMP_SCALE * tx_sams[i].real;

            fwrite(tx_scaled, sizeof(short), Nsamperframe, fout);
            frames++;
        }
        
	/* if this is in a pipeline, we probably don't want the usual
	   buffering to occur */

        if (fout == stdout) fflush(stdout);
        if (fin == stdin) fflush(stdin);

        if (testframes && (frames >= Nframes)) {
            goto finished;
        }
    }

 finished:
    fclose(fin);
    fclose(fout);

    ofdm_destroy(ofdm);

    return 0;
}

