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
#include "octave.h"
#include "mpdecode_core.h"
#include "gp_interleaver.h"
#include "interldpc.h"

#define NFRAMES  100               /* just log the first 100 frames          */
#define NDISCARD 20                /* BER2measure disctrds first 20 frames   */

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
    FILE          *fin, *fout, *foct;
    struct OFDM   *ofdm;
    int            nin_frame;

    float          phase_est_pilot_log[OFDM_ROWSPERFRAME*NFRAMES][OFDM_NC];
    COMP           rx_np_log[OFDM_ROWSPERFRAME*OFDM_NC*NFRAMES];
    float          rx_amp_log[OFDM_ROWSPERFRAME*OFDM_NC*NFRAMES];
    float          foff_hz_log[NFRAMES], snr_est_log[NFRAMES];
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
        fprintf(stderr, "  --testframes  Receive test frames and count errors\n");
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
    if (opt_exists(argv, argc, "--testframes")) {
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
    int iter[interleave_frames];
    int parityCheckCount[interleave_frames];
    
    for(i=0; i<interleave_frames; i++) {
        Nerrs_raw[i] = Nerrs_coded[i] = iter[i] = parityCheckCount[i] = 0;
    }
    
    ofdm = ofdm_create(OFDM_CONFIG_700D);
    assert(ofdm != NULL);

    if ((arg = opt_exists(argv, argc, "-v")) != 0) {
        ofdm_set_verbose(ofdm, 1);
    }
    if ((arg = opt_exists(argv, argc, "-vv")) != 0) {
        ofdm_set_verbose(ofdm, 2);
    }

    int Nbitsperframe = ofdm_get_bits_per_frame(ofdm);
    int Nmaxsamperframe = ofdm_get_max_samples_per_frame();
    
    short  rx_scaled[Nmaxsamperframe];
    COMP   rxbuf_in[Nmaxsamperframe];
    int    rx_bits[Nbitsperframe];
    char   rx_bits_char[Nbitsperframe];
    int    rx_uw[OFDM_NUWBITS];
    short  txt_bits[OFDM_NTXTBITS];
    f = 0; Nerrs = Terrs = Tbits = Terrs2 = Tbits2 = Terrs_coded = Tbits_coded = frame_count = 0;

    float EsNo = 3;
    fprintf(stderr,"Warning EsNo: %f hard coded\n", EsNo);

    float snr_est_smoothed_dB = 0.0;

    COMP  payload_syms[coded_syms_per_frame];
    float payload_amps[coded_syms_per_frame];
    COMP  codeword_symbols[interleave_frames*coded_syms_per_frame];
    float codeword_amps[interleave_frames*coded_syms_per_frame];
    for (i=0; i<interleave_frames*coded_syms_per_frame; i++) {
        codeword_symbols[i].real = 0.0;
        codeword_symbols[i].imag = 0.0;
        codeword_amps[i] = 0.0;
    }
    
    nin_frame = ofdm_get_nin(ofdm);
    while(fread(rx_scaled, sizeof(short), nin_frame, fin) == nin_frame) {

	/* scale and demod */

	for(i=0; i<nin_frame; i++) {
	    rxbuf_in[i].real = (float)rx_scaled[i]/(OFDM_AMP_SCALE/2);
            rxbuf_in[i].imag = 0.0;
        }

        if (strcmp(ofdm->sync_state,"search") == 0) {
            ofdm_sync_search(ofdm, rxbuf_in);
        }
        
        if ((strcmp(ofdm->sync_state,"synced") == 0) || (strcmp(ofdm->sync_state,"trial") == 0) ) {
            ofdm_demod(ofdm, rx_bits, rxbuf_in);
            ofdm_disassemble_modem_frame(ofdm, rx_uw, payload_syms, payload_amps, txt_bits);
            
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

                /* newest symbols at end of buffer (uses final i from last loop) */
                
                for(i=(interleave_frames-1)*coded_syms_per_frame,j=0; i<interleave_frames*coded_syms_per_frame; i++,j++) {
                    codeword_symbols[i] = payload_syms[j];
                    codeword_amps[i]    = payload_amps[j];
                }
               
                /* run de-interleaver */
                
                COMP  codeword_symbols_de[interleave_frames*coded_syms_per_frame];
                float codeword_amps_de[interleave_frames*coded_syms_per_frame];
                gp_deinterleave_comp (codeword_symbols_de, codeword_symbols, interleave_frames*coded_syms_per_frame);
                gp_deinterleave_float(codeword_amps_de   , codeword_amps   , interleave_frames*coded_syms_per_frame);

                double llr[coded_bits_per_frame];

                if (ldpc_en) {
                    char out_char[coded_bits_per_frame];

                    interleaver_sync_state_machine(ofdm, &ldpc, codeword_symbols_de, codeword_amps_de, EsNo,
                                                   interleave_frames, iter, parityCheckCount, Nerrs_coded);
                                         
                    if (!strcmp(ofdm->sync_state_interleaver,"synced") && (ofdm->frame_count_interleaver == interleave_frames)) {
                        ofdm->frame_count_interleaver = 0;
                        // printf("decode!\n");

                        if (testframes) {
                            Terrs += count_uncoded_errors(&ldpc, Nerrs_raw, interleave_frames, codeword_symbols_de);
                            Tbits += Nbitsperframe*interleave_frames;
                        }

                        for (j=0; j<interleave_frames; j++) {
                            symbols_to_llrs(llr, &codeword_symbols_de[j*coded_syms_per_frame],
                                                 &codeword_amps_de[j*coded_syms_per_frame],
                                                 EsNo, ofdm->mean_amp, coded_syms_per_frame);               
                            iter[j] = run_ldpc_decoder(&ldpc, out_char, llr, &parityCheckCount[j]);
                            //fprintf(stderr,"j: %d iter: %d pcc: %d\n", j, iter[j], parityCheckCount[j]);
                            if (testframes) {
                                Nerrs_coded[j] = count_errors(payload_data_bits, out_char, data_bits_per_frame);
                                Terrs_coded += Nerrs_coded[j];
                                Tbits_coded += data_bits_per_frame;
                            }
                            fwrite(out_char, sizeof(char), data_bits_per_frame, fout);
                        }
                    } /* if interleaver synced ..... */

                    /* SNR estimation and smoothing */

                    float snr_est_dB = 10*log10((ofdm->sig_var/ofdm->noise_var)*OFDM_NC*OFDM_RS/3000);
                    snr_est_smoothed_dB = 0.9*snr_est_smoothed_dB + 0.1*snr_est_dB;
                    
                } else {
                    /* lpdc_en == 0,  external LDPC decoder, so output LLRs */
                    symbols_to_llrs(llr, codeword_symbols_de, codeword_amps_de, EsNo, ofdm->mean_amp, coded_syms_per_frame);
                    fwrite(llr, sizeof(double), coded_bits_per_frame, fout);
                }
            } else {
                /* simple hard decision output for uncoded testing, all bits in frame dumped inlcuding UW and txt */
                for(i=0; i<Nbitsperframe; i++) {
                    rx_bits_char[i] = rx_bits[i];
                }
                fwrite(rx_bits_char, sizeof(char), Nbitsperframe, fout);
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
            int  r=0;
            if (testframes) {
                r = (ofdm->frame_count_interleaver - 1 ) % interleave_frames;
            }
            fprintf(stderr, "%3d st: %-6s euw: %2d %1d f: %5.1f ist: %-6s %2d eraw: %3d ecdd: %3d iter: %3d pcc: %3d",
                    f, ofdm->last_sync_state, ofdm->uw_errors, ofdm->sync_counter, ofdm->foff_est_hz,
                    ofdm->last_sync_state_interleaver, ofdm->frame_count_interleaver,
                    Nerrs_raw[r], Nerrs_coded[r], iter[r], parityCheckCount[r]);
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
            
            snr_est_log[f] = snr_est_smoothed_dB;
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
        octave_save_float(foct, "snr_est_log_c", snr_est_log, NFRAMES, 1, 1);
        fclose(foct);
    }

    ofdm_destroy(ofdm);

    return 0;
}
