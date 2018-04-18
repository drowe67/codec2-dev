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
  it under the terms of the GNU Lesser General Public License version 2, as
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
#include "test_bits_ofdm.h"

#define ASCALE (2E5*1.1491)

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
    int           testframes, ldpc_en, interleaver_frames;
   
    /* Set up default LPDC code.  We could add other codes here if we like */
    
    struct LDPC   ldpc;
    set_up_hra_112_112(&ldpc);
    int data_bits_per_frame = ldpc.data_bits_per_frame;
    int coded_bits_per_frame = ldpc.coded_bits_per_frame;
    int coded_syms_per_frame = ldpc.coded_syms_per_frame;
    
    if (argc < 3) {
        fprintf(stderr, "\n");
	fprintf(stderr, "usage: %s InputOneCharPerBitFile OutputModemRawFile [--lpdc] [--interleaver depth]\n\n", argv[0]);
        fprintf(stderr, "  -t            Transmit test frames (adjusts test frames for raw and LDPC modes)\n");
        fprintf(stderr, "  --ldpc        Run (%d,%d) LDPC decoder.  This forces 112, one char/bit output values\n"
                        "                per frame.  In testframe mode (-t) raw and coded errors will be counted\n",
                                         coded_bits_per_frame, data_bits_per_frame);
        fprintf(stderr, "  --interleave  Interleaver for LDPC frames, e.g. 1,2,4,8,16, default is 1\n");
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

    testframes = 0;
    if (opt_exists(argv, argc, "-t")) {
        testframes = 1;
    }

    /* set for LDPC coded or uncoded frames */
    
    ldpc_en = 0; interleaver_frames = 1;
    int Nbitsperframe;
    if (opt_exists(argv, argc, "--ldpc")) {

        assert((OFDM_NUWBITS+OFDM_NTXTBITS+coded_bits_per_frame) == OFDM_BITSPERFRAME); /* sanity check */

        ldpc_en = 1;
        if ((arg = opt_exists(argv, argc, "--interleaver"))) {
            interleaver_frames = atoi(argv[arg]);
        }
        Nbitsperframe = interleaver_frames*data_bits_per_frame;
        
    } else {
        /* vanilla uncoded input bits mode */
        Nbitsperframe = ofdm_get_bits_per_frame(ofdm);
    }
    
    int Nsamperframe = ofdm_get_samples_per_frame();

    unsigned char tx_bits_char[Nbitsperframe];
    int           tx_bits[Nbitsperframe];
    short         tx_scaled[Nsamperframe];
    
    /* build modulated UW and txt bits */

    int  uw_txt_bits[OFDM_NUWBITS+OFDM_NTXTBITS];
    COMP uw_txt_syms[(OFDM_NUWBITS+OFDM_NTXTBITS)/OFDM_BPS];
    complex float tx_symbols[(OFDM_NUWBITS+OFDM_NTXTBITS)/OFDM_BPS + coded_syms_per_frame];
    for(i=0; i<OFDM_NUWBITS; i++) {
        uw_txt_bits[i] = ofdm->tx_uw[i];
    }
    for(j=0; j<OFDM_NTXTBITS; j++,i++) {
        uw_txt_bits[i] = 0;
    }    
    qpsk_modulate_frame(uw_txt_syms, uw_txt_bits, (OFDM_NUWBITS+OFDM_NTXTBITS)/OFDM_BPS);
    for(i=0; i<(OFDM_NUWBITS+OFDM_NTXTBITS)/OFDM_BPS; i++) {
        tx_symbols[i] = uw_txt_syms[i].real + I * uw_txt_syms[i].imag;
    }
    
    /* main loop ----------------------------------------------------------------*/
    
    frames = 0;

    while(fread(tx_bits_char, sizeof(char), Nbitsperframe, fin) == Nbitsperframe) {
	frames++;

        if (ldpc_en) {
            /* fancy interleaved LDPC encoded frames ----------------------------------------*/
            
            /* optionally overwrite input data with test frame nown to demodulator */
            
            if (testframes) {
                for (j=0; j<interleaver_frames; j++) {
                    for(i=0; i<data_bits_per_frame; i++) {
                        tx_bits_char[j*data_bits_per_frame + i] = payload_data_bits[i];
                    }
                }
            }
            
            int codeword[coded_bits_per_frame];
            COMP coded_symbols[interleaver_frames*coded_syms_per_frame];
            COMP coded_symbols_inter[interleaver_frames*coded_syms_per_frame];
            complex float tx_sams[Nsamperframe];

            for (j=0; j<interleaver_frames; j++) {
                ldpc_encode_frame(&ldpc, codeword, &tx_bits_char[j*data_bits_per_frame]);
                qpsk_modulate_frame(&coded_symbols[j*coded_syms_per_frame], codeword, coded_syms_per_frame);
                gp_interleave_comp(coded_symbols_inter, coded_symbols, interleaver_frames*coded_syms_per_frame);
                for(i=0; i<coded_syms_per_frame; i++) {
                    tx_symbols[(OFDM_NUWBITS+OFDM_NTXTBITS)/OFDM_BPS+i] = coded_symbols_inter[j*coded_syms_per_frame+i].real
                                                                        + I * coded_symbols_inter[j*coded_syms_per_frame+i].imag;
                }
                ofdm_txframe(ofdm, tx_sams, tx_symbols);
                for(i=0; i<Nsamperframe; i++) {
                    tx_scaled[i] = ASCALE * crealf(tx_sams[i]);
                }
                fwrite(tx_scaled, sizeof(short), Nsamperframe, fout);
            }
        
        } else {
            /* just modulate uncoded raw bits ----------------------------------------------*/
            
            if (testframes) {
                for(i=0; i<Nbitsperframe; i++) {
                    tx_bits_char[i] = test_bits_ofdm[i];
                }
            }
            for(i=0; i<Nbitsperframe; i++)
                tx_bits[i] = tx_bits_char[i];
            COMP tx_sams[Nsamperframe];
            ofdm_mod(ofdm, tx_sams, tx_bits);

            /* scale and save to disk as shorts */

            for(i=0; i<Nsamperframe; i++)
                tx_scaled[i] = ASCALE * tx_sams[i].real;

            fwrite(tx_scaled, sizeof(short), Nsamperframe, fout);
        }
        
	/* if this is in a pipeline, we probably don't want the usual
	   buffering to occur */

        if (fout == stdout) fflush(stdout);
        if (fin == stdin) fflush(stdin);
    }

    fclose(fin);
    fclose(fout);
    ofdm_destroy(ofdm);

    return 0;
}
