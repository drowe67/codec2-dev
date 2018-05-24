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
    int           frames;
    int           i, j, arg;
    int           testframes, ldpc_en, interleave_frames;
   
    /* Set up default LPDC code.  We could add other codes here if we like */
    
    struct LDPC   ldpc;
    set_up_hra_112_112(&ldpc);
    int data_bits_per_frame = ldpc.data_bits_per_frame;
    int coded_bits_per_frame = ldpc.coded_bits_per_frame;
    
    if (argc < 3) {
        fprintf(stderr, "\n");
	fprintf(stderr, "usage: %s InputOneCharPerBitFile OutputModemRawFile [--lpdc] [--interleaver depth]\n\n", argv[0]);
        fprintf(stderr, "  --testframes Nsecs  Transmit test frames (adjusts test frames for raw and LDPC modes)\n");
        fprintf(stderr, "  --ldpc              Run (%d,%d) LDPC decoder.  This forces 112, one char/bit output values\n"
                        "                      per frame.  In testframe mode (-t) raw and coded errors will be counted\n",
                                               coded_bits_per_frame, data_bits_per_frame);
        fprintf(stderr, "  --interleave depth  Interleave depth for LDPC frames, e.g. 1,2,4,8,16, default is 1\n");
        fprintf(stderr, "  --txbpf             Transmit band pass filter\n");
        fprintf(stderr, "\n");
	exit(1);
    }

    if (strcmp(argv[1], "-")  == 0) fin = stdin;
    else if ( (fin = fopen(argv[1],"rb")) == NULL ) {
        fprintf(stderr, "Error opening input file: %s: %s.\n",
                argv[1], strerror(errno));
        exit(1);
    }

    if (strcmp(argv[2], "-") == 0) fout = stdout;
    else if ( (fout = fopen(argv[2],"wb")) == NULL ) {
        fprintf(stderr, "Error opening output modem sample file: %s: %s.\n",
                argv[2], strerror(errno));
        exit(1);
    }
    
    ofdm = ofdm_create(OFDM_CONFIG_700D);
    assert(ofdm != NULL);

    /* set for LDPC coded or uncoded frames */
    
    ldpc_en = 0; interleave_frames = 1;
    int Nbitsperframe;
    if (opt_exists(argv, argc, "--ldpc")) {

        assert((OFDM_NUWBITS+OFDM_NTXTBITS+coded_bits_per_frame) == OFDM_BITSPERFRAME); /* sanity check */

        ldpc_en = 1;
        if ((arg = opt_exists(argv, argc, "--interleave"))) {
            interleave_frames = atoi(argv[arg+1]);
        }
        Nbitsperframe = interleave_frames*data_bits_per_frame;
        
    } else {
        /* vanilla uncoded input bits mode */
        Nbitsperframe = ofdm_get_bits_per_frame(ofdm);
    }
    
    int Nsamperframe = ofdm_get_samples_per_frame();
    fprintf(stderr, "Nbitsperframe: %d interleave_frames: %d\n", Nbitsperframe, interleave_frames);

    uint8_t tx_bits_char[Nbitsperframe];
    short   tx_scaled[Nsamperframe];
    uint8_t txt_bits_char[OFDM_NTXTBITS*interleave_frames];

    for(i=0; i<OFDM_NTXTBITS*interleave_frames; i++) { txt_bits_char[i] = 0; }
   
    testframes = 0;
    int Nframes = 0;
    if ((arg = (opt_exists(argv, argc, "--testframes")))) {
        testframes = 1;
        int Nsec, Nrows;
        Nsec = atoi(argv[arg+1]);
        Nrows = Nsec*OFDM_RS;
        Nframes = floor((Nrows-1)/OFDM_NS);
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
                for (j=0; j<interleave_frames; j++) {
                    for(i=0; i<data_bits_per_frame; i++) {
                        tx_bits_char[j*data_bits_per_frame + i] = payload_data_bits[i];
                    }
                }
            }

            complex float tx_sams[interleave_frames*Nsamperframe];
            ofdm_ldpc_interleave_tx(ofdm, &ldpc, tx_sams, tx_bits_char, txt_bits_char, interleave_frames);

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
                for(i=0; i<Nbitsperframe; i++) {
                    tx_bits_char[i] = test_bits_ofdm[i];
                }
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
