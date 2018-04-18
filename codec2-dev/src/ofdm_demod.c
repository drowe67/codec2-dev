/*---------------------------------------------------------------------------*\

  FILE........: ofdm_demod.c
  AUTHOR......: David Rowe
  DATE CREATED: Mar 2018

  Given an input file of raw file (8kHz, 16 bit shorts) of OFDM modem
  samples.  Optionally:

    1/ outputs one char per bit (hard decision)
    2/ bit LLRS, one double per bir, for external LDPC decoder like ldpc_dec
    3/ LDPC decoded bits, one char per bit

  Also has test frame modes for uncoded and coded operation.

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
#include "octave.h"
#include "test_bits_ofdm.h"
#include "mpdecode_core.h"
#include "gp_interleaver.h"
#include "interldpc.h"

#define ASCALE   (2E5*1.1491/2.0)  /* scale from shorts back to floats       */
#define NFRAMES  100               /* just log the first 100 frames          */
#define NDISCARD 20                /* BER2measure disctrds first 20 frames   */

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
    FILE          *fin, *fout, *foct;
    struct OFDM   *ofdm;
    int            nin_frame;

    float          phase_est_pilot_log[OFDM_ROWSPERFRAME*NFRAMES][OFDM_NC];
    COMP           rx_np_log[OFDM_ROWSPERFRAME*OFDM_NC*NFRAMES];
    float          rx_amp_log[OFDM_ROWSPERFRAME*OFDM_NC*NFRAMES];
    float          foff_hz_log[NFRAMES];
    int            timing_est_log[NFRAMES];

    int            i, j, f, oct, logframes, arg, llr_en, interleave_frames;
    int            Nerrs, Terrs, Tbits, Terrs2, Tbits2, testframes, frame_count;
    int            ldpc_en, Tbits_coded, Terrs_coded;
    
    /* Set up default LPDC code.  We could add other codes here if we like */
    
    struct LDPC ldpc;
    set_up_hra_112_112(&ldpc);
    int data_bits_per_frame = ldpc.data_bits_per_frame;
    int coded_bits_per_frame = ldpc.coded_bits_per_frame;
    int coded_syms_per_frame = ldpc.coded_syms_per_frame;

    if (argc < 3) {
        fprintf(stderr, "\n");
	printf("usage: %s InputModemRawFile OutputFile [-o OctaveLogFile] [--llr] [--ldpc] [--interleave depth] [-v]\n", argv[0]);
        fprintf(stderr, "\n");
        fprintf(stderr, "                Default output file format is one byte per bit hard decision\n");
        fprintf(stderr, "  --llr         LLR output, one double per bit, %d doubles/frame\n", coded_bits_per_frame);
        fprintf(stderr, "  -t            Receive test frames and count errors\n");
        fprintf(stderr, "  --ldpc        Run (%d,%d) LDPC decoder.  This forces 112, one char/bit output values\n"
                        "                per frame.  In testframe mode (-t) raw and coded errors will be counted\n",
                                         coded_bits_per_frame, data_bits_per_frame);
        fprintf(stderr, "  --interleave  Interleaver for LDPC frames, e.g. 1,2,4,8,16, default is 1\n");
        fprintf(stderr, "  -v            Verbose info the stderr\n");
        fprintf(stderr, "  -o            Octave log file for testing\n");
        fprintf(stderr, "\n");
	exit(1);
    }

    if (strcmp(argv[1], "-")  == 0) fin = stdin;
    else if ( (fin = fopen(argv[1],"rb")) == NULL ) {
	fprintf(stderr, "Error opening input modem sample file: %s: %s.\n",
         argv[1], strerror(errno));
	exit(1);
    }

    if (strcmp(argv[2], "-") == 0) fout = stdout;
    else if ( (fout = fopen(argv[2],"wb")) == NULL ) {
	fprintf(stderr, "Error opening output file: %s: %s.\n",
         argv[2], strerror(errno));
	exit(1);
    }

    foct = NULL;
    oct = 0;
    if ((arg = opt_exists(argv, argc, "-o")) != 0) {
        if ( (foct = fopen(argv[arg+1],"wt")) == NULL ) {
            fprintf(stderr, "Error opening output Octave file: %s: %s.\n",
                    argv[4], strerror(errno));
	exit(1);
        }
        oct = 1;
        logframes = NFRAMES;
    }

    llr_en = 0;
    if (opt_exists(argv, argc, "--llr")) {
        llr_en = 1;
    }

    testframes = 0;
    if (opt_exists(argv, argc, "-t")) {
        testframes = 1;
    }

    ldpc_en = 0;
    if (opt_exists(argv, argc, "--ldpc")) {
        ldpc_en = 1;
        llr_en = 1;
    }

    interleave_frames = 1;
    if ((arg = opt_exists(argv, argc, "--interleave"))) {
        interleave_frames = atoi(argv[arg+1]);
        fprintf(stderr, "interleave_frames: %d\n",  interleave_frames);
        /* we can't de-interleave without LDPC for sync, so switch that on */
        ldpc_en = 1;
        llr_en = 1;
    }
    int Nerrs_raw[interleave_frames];
    int Nerrs_coded[interleave_frames];
    for(i=0; i<interleave_frames; i++) {
        Nerrs_raw[i] = Nerrs_coded[i] = 0;
    }
    
    ofdm = ofdm_create(OFDM_CONFIG_700D);
    assert(ofdm != NULL);

    if ((arg = opt_exists(argv, argc, "-v")) != 0) {
        ofdm_set_verbose(ofdm, 1);
    }

    int Nbitsperframe = ofdm_get_bits_per_frame(ofdm);
    int Nmaxsamperframe = ofdm_get_max_samples_per_frame();
    
    short  rx_scaled[Nmaxsamperframe];
    COMP   rxbuf_in[Nmaxsamperframe];
    int    rx_bits[Nbitsperframe];
    char   rx_bits_char[Nbitsperframe];
    int    rx_uw[OFDM_NUWBITS];
    f = 0; Nerrs = Terrs = Tbits = Terrs2 = Tbits2 = Terrs_coded = Tbits_coded = frame_count = 0;
    int    parityCheckCount, iter;

    float EsNo = 10;
    fprintf(stderr,"Warning EsNo: %f hard coded\n", EsNo);

    COMP  codeword_symbols[interleave_frames*coded_syms_per_frame];
    float codeword_amps[interleave_frames*coded_syms_per_frame];
    COMP  codeword_symbols_de[interleave_frames*coded_syms_per_frame];
    float codeword_amps_de[interleave_frames*coded_syms_per_frame];

    nin_frame = ofdm_get_nin(ofdm);
    while(fread(rx_scaled, sizeof(short), nin_frame, fin) == nin_frame) {

	/* scale and demod */

	for(i=0; i<nin_frame; i++) {
	    rxbuf_in[i].real = (float)rx_scaled[i]/ASCALE;
            rxbuf_in[i].imag = 0.0;
        }

        if (strcmp(ofdm->sync_state,"search") == 0) {
            ofdm_sync_search(ofdm, rxbuf_in);
        }
    
        if ((strcmp(ofdm->sync_state,"synced") == 0) || (strcmp(ofdm->sync_state,"trial") == 0) ) {
            ofdm_demod(ofdm, rx_bits, rxbuf_in);
            
            if (llr_en) {
                
                /* first few symbols are used for UW and txt bits, find start of (224,112) LDPC codeword 
                   and extract QPSK symbols and amplitude estimates */

                assert((OFDM_NUWBITS+OFDM_NTXTBITS+coded_bits_per_frame) == OFDM_BITSPERFRAME);

                /* now we need to buffer for de-interleaving -------------------------------------*/
                
                /* shift interleaved symbol buffers to make room for new symbols */
                
                for(i=0, j=coded_syms_per_frame; j<interleave_frames*coded_syms_per_frame; i++,j++) {
                    codeword_symbols[i] = codeword_symbols[j];
                    codeword_amps[i] = codeword_amps[j];
                }

                /* newest symbols at end of buffer (uses final i from last loop), note we 
                   change COMP formats from what modem uses internally */
                
                for(i=(interleave_frames-1)*coded_syms_per_frame,j=(OFDM_NUWBITS+OFDM_NTXTBITS)/OFDM_BPS; i<interleave_frames*coded_syms_per_frame; i++,j++) {
                    codeword_symbols[i].real = crealf(ofdm->rx_np[j]);
                    codeword_symbols[i].imag = cimagf(ofdm->rx_np[j]);
                    codeword_amps[i] = ofdm->rx_amp[j];
                }
               
                /* run de-interleaver */
                
                gp_deinterleave_comp (codeword_symbols_de, codeword_symbols, interleave_frames*coded_syms_per_frame);
                gp_deinterleave_float(codeword_amps_de   , codeword_amps   , interleave_frames*coded_syms_per_frame);

                double llr[coded_bits_per_frame];

                if (ldpc_en) {
                    char out_char[coded_bits_per_frame];

                    /* 
                       Interleaver Sync:
                         Needs to work on any data
                         Use indication of LDPC convergence, may need to patch CML code for that
                         Attempt a decode on every frame, when it converges we have sync
                    */
                    
                    char next_sync_state_interleaver[OFDM_STATE_STR];
                    strcpy(next_sync_state_interleaver, ofdm->sync_state_interleaver);
                    if (strcmp(ofdm->sync_state_interleaver,"search") == 0) {
                        symbols_to_llrs(llr, codeword_symbols_de, codeword_amps_de, EsNo, coded_syms_per_frame);               
                        iter = run_ldpc_decoder(&ldpc, out_char, llr, &parityCheckCount);
                        Nerrs = data_bits_per_frame - parityCheckCount;
                        //fprintf(stderr, "iter: %d pcc: %d Nerrs: %d\n", iter, parityCheckCount, Nerrs);
                        if (Nerrs < 10) {
                            /* sucessful decode! */
                            strcpy(next_sync_state_interleaver, "synced");
                            ofdm->frame_count_interleaver = interleave_frames;
                        }
                    }
                    strcpy(ofdm->sync_state_interleaver, next_sync_state_interleaver);
                     
                    if (!strcmp(ofdm->sync_state_interleaver,"synced") && (ofdm->frame_count_interleaver == interleave_frames)) {
                        ofdm->frame_count_interleaver = 0;
                        // printf("decode!\n");

                        if (testframes) {
                            
                            /* measure uncoded (raw) bit errors over interleaver frame */

                            int rx_bits_raw[coded_bits_per_frame];
                            for (j=0; j<interleave_frames; j++) {
                                for(i=0; i<coded_syms_per_frame; i++) {
                                    int bits[2];
                                    complex float s = codeword_symbols_de[j*coded_syms_per_frame+i].real + I*codeword_symbols_de[j*coded_syms_per_frame+i].imag;
                                    qpsk_demod(s, bits);
                                    rx_bits_raw[OFDM_BPS*i]   = bits[1];
                                    rx_bits_raw[OFDM_BPS*i+1] = bits[0];
                                }
                                Nerrs = 0;
                                assert(sizeof(test_codeword)/sizeof(int) == coded_bits_per_frame);
                                for(i=0; i<coded_bits_per_frame; i++) {
                                    //fprintf(stderr, "%d %d %d\n", i, test_codeword[i], rx_bits_raw[i]);
                                    if (test_codeword[i] != rx_bits_raw[i]) {
                                        Nerrs++;
                                    }
                                }
                                
                                Nerrs_raw[j] = Nerrs;
                                Terrs += Nerrs;
                                Tbits += Nbitsperframe;
                            }
                        }

                        for (j=0; j<interleave_frames; j++) {
                            symbols_to_llrs(llr, &codeword_symbols_de[j*coded_syms_per_frame],
                                                 &codeword_amps_de[j*coded_syms_per_frame],
                                                 EsNo, coded_syms_per_frame);               
                            iter = run_ldpc_decoder(&ldpc, out_char, llr, &parityCheckCount);

                            if (testframes) {
                                Nerrs = 0;
                                for(i=0; i<data_bits_per_frame; i++) {
                                    if (payload_data_bits[i] != out_char[i]) {
                                        Nerrs++;
                                    }
                                }
                                Nerrs_coded[j] = Nerrs;
                                Terrs_coded += Nerrs;
                                Tbits_coded += data_bits_per_frame;
                            }
                            fwrite(out_char, sizeof(char), coded_bits_per_frame, fout);
                        }
                    } /* if interleaver synced ..... */
                         
                } else {
                    /* lpdc_en == 0,  external LDPC decoder, so output LLRs */
                    symbols_to_llrs(llr, codeword_symbols_de, codeword_amps_de, EsNo, coded_syms_per_frame);
                        fwrite(llr, sizeof(double), coded_bits_per_frame, fout);
                }
            } else {
                /* simple hard decision output for uncoded testing, all bits in frame dumped inlcuding UW and txt */
                for(i=0; i<Nbitsperframe; i++) {
                    rx_bits_char[i] = rx_bits[i];
                }
                fwrite(rx_bits_char, sizeof(char), Nbitsperframe, fout);
            }

            /* extract Unique Word bits */

            for(i=0; i<OFDM_NUWBITS; i++) {
                rx_uw[i] = rx_bits[i];
            }

            /* optional error counting on uncoded data in non-LDPC testframe mode */
            
            if (testframes && (ldpc_en == 0)) {
                Nerrs = 0;
                for(i=0; i<Nbitsperframe; i++) {
                    if (test_bits_ofdm[i] != rx_bits[i]) {
                        Nerrs++;
                    }
                }

                Terrs += Nerrs;
                Tbits += Nbitsperframe;

                if (frame_count >= NDISCARD) {
                    Terrs2 += Nerrs;
                    Tbits2 += Nbitsperframe;
                }
            }
            frame_count++;
        }
        
        nin_frame = ofdm_get_nin(ofdm);
        ofdm_sync_state_machine(ofdm, rx_uw);

        /* act on any events returned by state machine */
    
        if (ofdm->sync_start) {
            Terrs = Tbits = Terrs2 = Tbits2 = Terrs_coded = Tbits_coded = frame_count = 0;
            for(i=0; i<interleave_frames; i++) {
                Nerrs_raw[i] = Nerrs_coded[i] = 0;
            }

        }

        if (ofdm->verbose) {
           int  r = ofdm->frame_count_interleaver % interleave_frames;
            fprintf(stderr, "f: %3d st: %-6s uw_errs: %2d %1d inter_st: %-6s inter_fr: %2d Nerrs_raw: %3d Nerrs_coded: %3d foff: %4.1f",
                    f, ofdm->last_sync_state, ofdm->uw_errors, ofdm->sync_counter,
                    ofdm->last_sync_state_interleaver, ofdm->frame_count_interleaver,
                    Nerrs_raw[r], Nerrs_coded[r], ofdm->foff_est_hz);
            fprintf(stderr, "\n");    
        }

        /* optional logging of states */
        
        if (oct) {
            /* note corrected phase (rx no phase) is one big linear array for frame */

            for (i = 0; i < OFDM_ROWSPERFRAME*OFDM_NC; i++) {
                rx_np_log[OFDM_ROWSPERFRAME*OFDM_NC*f + i].real = crealf(ofdm->rx_np[i]);
                rx_np_log[OFDM_ROWSPERFRAME*OFDM_NC*f + i].imag = cimagf(ofdm->rx_np[i]);
            }

            /* note phase/amp ests the same for each col, but check them all anyway */

            for (i = 0; i < OFDM_ROWSPERFRAME; i++) {
                for (j = 0; j < OFDM_NC; j++) {
                    phase_est_pilot_log[OFDM_ROWSPERFRAME*f+i][j] = ofdm->aphase_est_pilot_log[OFDM_NC*i+j];
                    rx_amp_log[OFDM_ROWSPERFRAME*OFDM_NC*f+OFDM_NC*i+j] = ofdm->rx_amp[OFDM_NC*i+j];
                }
            }

            foff_hz_log[f] = ofdm->foff_est_hz;
            timing_est_log[f] = ofdm->timing_est + 1;     /* offset by 1 to match Octave */
            if (f == (logframes-1))
                oct = 0;
        }

	/* if this is in a pipeline, we probably don't want the usual
	   buffering to occur */

        if (fout == stdout) fflush(stdout);
        if (fin == stdin) fflush(stdin);

        f++;
    }

    fclose(fin);
    fclose(fout);

    if (testframes) {
        fprintf(stderr, "BER......: %5.4f Tbits: %5d Terrs: %5d\n", (float)Terrs/Tbits, Tbits, Terrs);
        if (!ldpc_en) {
            fprintf(stderr, "BER2.....: %5.4f Tbits: %5d Terrs: %5d\n", (float)Terrs2/Tbits2, Tbits2, Terrs2);
        }
        if (ldpc_en) {
            fprintf(stderr, "Coded BER: %5.4f Tbits: %5d Terrs: %5d\n",
                    (float)Terrs_coded/Tbits_coded, Tbits_coded, Terrs_coded);
        }        
    }
    
    /* optionally dump Octave files */

    if (foct != NULL) {
        octave_save_float(foct, "phase_est_pilot_log_c", (float*)phase_est_pilot_log, OFDM_ROWSPERFRAME*NFRAMES, OFDM_NC, OFDM_NC);
        octave_save_complex(foct, "rx_np_log_c", (COMP*)rx_np_log, 1, OFDM_ROWSPERFRAME*OFDM_NC*NFRAMES, OFDM_ROWSPERFRAME*OFDM_NC*NFRAMES);
        octave_save_float(foct, "rx_amp_log_c", (float*)rx_amp_log, 1, OFDM_ROWSPERFRAME*OFDM_NC*NFRAMES, OFDM_ROWSPERFRAME*OFDM_NC*NFRAMES);
        octave_save_float(foct, "foff_hz_log_c", foff_hz_log, NFRAMES, 1, 1);
        octave_save_int(foct, "timing_est_log_c", timing_est_log, NFRAMES, 1);
        fclose(foct);
    }

    ofdm_destroy(ofdm);

    return 0;
}
