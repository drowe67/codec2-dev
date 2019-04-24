/*---------------------------------------------------------------------------*\

  FILE........: tst_ofdm_demod.c
  AUTHOR......: David Rowe, Don Reid
  DATE CREATED: 7 July 2018

  Test and profile OFDM de-modulation on the STM32F4.

\*---------------------------------------------------------------------------*/

/*
  Copyright (C) 2014 David Rowe

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


/* This is a unit test implementation of the OFDM Demod function.
 * It is used for several tests:
 *
 *  tst_ofdm_demod_ideal    Simple 10 frames with no degradation.
 *  tst_ofdm_demod_AWGN     Just AWGN in channel.
 *  tst_ofdm_demod_fade     AWGN and fading in channel.
 *  tst_ofdm_demod_profile  Profile, disable verbose logging.
 *
 * See tst_ofdm_demod_setup and tst_ofdm_demod_check scripts for details.
 *
 * This program reads a file "stm_cfg.txt" at startup to configure its options.
 *
 * This program is intended to be run using input data, typically
 * Codec2 frames, which may have had simulated RF degredation applied.
 * For example:
 *
 *    ofdm_get_test_bits - 10 | * ofdm_mod - - | \
 *        cohpsk_ch - stm_in.raw -20 -Fs 8000 -f -5
 *
 * Reference data can be created by running the same input through the x86 
 * ofdm_demod tool.
 *
 *    ofdm_demod stm_in.raw ref_demod_out.raw -o ofdm_demod_ref_log.txt --testframes
 *
 * Comparison of the results to the reference will depend on the test conditions.
 * Some small differences are expected due to differences in implementation.
 *
 */


#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <errno.h>
#include <fcntl.h>
#include <unistd.h>

#include "semihosting.h"
#include "codec2_ofdm.h"
#include "ofdm_internal.h"
#include "mpdecode_core.h"
#include "interldpc.h"
#include "gp_interleaver.h"
#include "test_bits_ofdm.h"

#include "debug_alloc.h"

#include "stm32f4xx_conf.h"
#include "stm32f4xx.h"
#include "machdep.h"

#define NDISCARD 20

extern const uint8_t payload_data_bits[];
extern const int test_bits_ofdm[];

static struct OFDM_CONFIG *ofdm_config;

static int ofdm_bitsperframe;
static int ofdm_rowsperframe;
static int ofdm_nuwbits;
static int ofdm_ntxtbits;
static int ofdm_nin;
static char fout_buffer[4*4096];
static __attribute__ ((section (".ccm"))) char fdiag_buffer[4*8192];
static __attribute__ ((section (".ccm"))) char fin_buffer[4096*8];

static char *statemode[] = {
    "search",
    "trial",
    "synced"
};

static FILE *fout, *fdiag;
void flush_all(void) {
    fflush(fout);
    fflush(fdiag);
    fflush(stdout);
    fflush(stderr);
    }

int main(int argc, char *argv[]) {
    struct OFDM *ofdm;
    FILE        *fcfg;
    int          nin_frame;
    struct LDPC  ldpc;

    // Test configuration, read from stm_cfg.txt
    int          config_verbose; 
    int          config_testframes; 
    int          config_ldpc_en; 
    int          config_log_payload_syms;
    int          config_profile;

    int          i, j, f;
    int          interleave_frames;
    int          Nerrs, Terrs, Tbits, Terrs2, Tbits2, frame_count;
    int          Tbits_coded, Terrs_coded;

    semihosting_init();

    fprintf(stdout, "OFDM Demod test\n");

    // Read configuration - a file of '0' or '1' characters
    char config[8];
    fcfg = fopen("stm_cfg.txt", "r");
    if (fcfg == NULL) {
        fprintf(stderr, "Error opening config file\n");
        exit(1);
    }
    if (fread(&config[0], 1, 8, fcfg) != 8) {
        fprintf(stderr, "Error reading config file\n");
        exit(1);
    }
    config_verbose = config[0] - '0';
    config_testframes = config[1] - '0';
    config_ldpc_en = config[2] - '0';
    config_log_payload_syms = config[3] - '0';
    config_profile = config[4] - '0';
    fclose(fcfg);
    //
    interleave_frames = 1;

    int Nerrs_raw[interleave_frames];
    int Nerrs_coded[interleave_frames];
    int iter[interleave_frames];
    int parityCheckCount[interleave_frames];

    for(i=0; i<interleave_frames; i++) {
        Nerrs_raw[i] = Nerrs_coded[i] = iter[i] = parityCheckCount[i] = 0;
    }

    PROFILE_VAR(ofdm_demod_start, ofdm_demod_sync_search,
                ofdm_demod_demod, ofdm_demod_diss, ofdm_demod_snr);
    ofdm_demod_start = 0;
    ofdm_demod_sync_search = 0;
    ofdm_demod_demod = 0;
    ofdm_demod_diss = 0;
    ofdm_demod_snr = 0;
    if (config_profile) machdep_profile_init();

    if ((ofdm_config = (struct OFDM_CONFIG *) CALLOC(1, sizeof (struct OFDM_CONFIG))) == NULL) {
        fprintf(stderr, "Out of Memory\n");
        exit(1);
    }

    ofdm_config->fs = 8000.0;			/* Sample Frequency */
    ofdm_config->ofdm_timing_mx_thresh = 0.30;
    ofdm_config->ftwindowwidth = 11;
    ofdm_config->bps = 2;   			/* Bits per Symbol */
    ofdm_config->txtbits = 4; 			/* number of auxiliary data bits */
    ofdm_config->ns = 8;  			/* Number of Symbol frames */

    ofdm_config->rs = (1.0f / ofdm_config->ts); /* Symbol Rate */

    /* config options can change here, non yet */

    ofdm = ofdm_create(ofdm_config);
    assert(ofdm != NULL);

    FREE(ofdm_config);

    /* Get a copy of the actual modem config */
    ofdm_config = ofdm_get_config_param();

    set_up_hra_112_112(&ldpc, ofdm_config);

    ofdm_bitsperframe = ofdm_get_bits_per_frame();
    ofdm_rowsperframe = ofdm_bitsperframe / (ofdm_config->nc * ofdm_config->bps);
    ofdm_nuwbits = (ofdm_config->ns - 1) * ofdm_config->bps - ofdm_config->txtbits;
    ofdm_ntxtbits = ofdm_config->txtbits;
    ofdm_nin = ofdm_get_nin(ofdm);

    ofdm_set_verbose(ofdm, config_verbose);

    int Nmaxsamperframe = ofdm_get_max_samples_per_frame();

    int data_bits_per_frame = ldpc.data_bits_per_frame;
    int coded_bits_per_frame = ldpc.coded_bits_per_frame;
    int coded_syms_per_frame = ldpc.coded_syms_per_frame;

    short   rx_scaled[Nmaxsamperframe];
    int     rx_bits[ofdm_bitsperframe];
    char    rx_bits_char[ofdm_bitsperframe];
    uint8_t rx_uw[ofdm_nuwbits];
    short   txt_bits[ofdm_ntxtbits];
    f = 0;
    Nerrs = Terrs = Tbits = Terrs2 = Tbits2 = Terrs_coded = Tbits_coded = frame_count = 0;

    float snr_est_smoothed_dB = 0.0;

    float EsNo = 3.0f;  // Constant from ofdm_demod.c

    COMP  payload_syms[coded_syms_per_frame];
    float payload_amps[coded_syms_per_frame];
    COMP  codeword_symbols[interleave_frames*coded_syms_per_frame];
    float codeword_amps[interleave_frames*coded_syms_per_frame];

    FILE* fin = fopen("stm_in.raw", "rb");
    if (fin == NULL) {
        fprintf(stderr, "Error opening input file\n");
        exit(1);
    }
    setvbuf(fin, fin_buffer,_IOFBF,sizeof(fin_buffer));


    fout = fopen("stm_out.raw", "wb");
    if (fout == NULL) {
        fprintf(stderr, "Error opening output file\n");
        exit(1);
    }
    setvbuf(fout, fout_buffer,_IOFBF,sizeof(fout_buffer));

    fdiag = fopen("stm_diag.raw", "wb");
    if (fdiag == NULL) {
        fprintf(stderr, "Error opening diag file\n");
        exit(1);
    }
    setvbuf(fdiag, fdiag_buffer,_IOFBF,sizeof(fdiag_buffer));

    nin_frame = ofdm_get_nin(ofdm);
    int num_read;

    while((num_read = fread(rx_scaled, sizeof(short) , nin_frame, fin)) == nin_frame) {

        int log_payload_syms_flag = 0;

        if (config_profile) PROFILE_SAMPLE(ofdm_demod_start);

	    /* demod */

        if (config_profile) PROFILE_SAMPLE_AND_LOG2(ofdm_demod_start, "  ofdm_demod_start");

        if (ofdm->sync_state == search) {
            if (config_profile) PROFILE_SAMPLE(ofdm_demod_sync_search);
            ofdm_sync_search_shorts(ofdm, rx_scaled, (OFDM_AMP_SCALE/2));
            if (config_profile) PROFILE_SAMPLE_AND_LOG2(ofdm_demod_sync_search, "  ofdm_demod_sync_search");
        }

        if ((ofdm->sync_state == synced) || (ofdm->sync_state == trial) ) {
            if (config_profile) PROFILE_SAMPLE(ofdm_demod_demod);
            ofdm_demod_shorts(ofdm, rx_bits, rx_scaled, (OFDM_AMP_SCALE/2));
            if (config_profile) PROFILE_SAMPLE_AND_LOG2(ofdm_demod_demod, "  ofdm_demod_demod");
            if (config_profile) PROFILE_SAMPLE(ofdm_demod_diss);
            ofdm_disassemble_modem_frame(ofdm, rx_uw, payload_syms, payload_amps, txt_bits);
            if (config_profile) PROFILE_SAMPLE_AND_LOG2(ofdm_demod_diss, "  ofdm_demod_diss");
            log_payload_syms_flag = 1;

            /* SNR estimation and smoothing */
            if (config_profile) PROFILE_SAMPLE(ofdm_demod_snr);
            float snr_est_dB = 10*log10((ofdm->sig_var/ofdm->noise_var) * 
                                ofdm_config->nc * ofdm_config->rs / 3000);
            snr_est_smoothed_dB = 0.9f * snr_est_smoothed_dB + 0.1f  *snr_est_dB;
            if (config_profile) {
                PROFILE_SAMPLE_AND_LOG2(ofdm_demod_snr, "  ofdm_demod_snr");
            }

            // LDPC
            if (config_ldpc_en) {  // was llr_en in orig

                /* first few symbols are used for UW and txt bits, find
                   start of (224,112) LDPC codeword and extract QPSK
                   symbols and amplitude estimates */
                assert((ofdm_nuwbits + ofdm_ntxtbits + coded_bits_per_frame) 
                        == ofdm_bitsperframe);

                /* now we need to buffer for de-interleaving */

                /* shift interleaved symbol buffers to make room for new symbols */
                for(i=0, j=coded_syms_per_frame; 
                    j<interleave_frames*coded_syms_per_frame; 
                    i++,j++) {
                    codeword_symbols[i] = codeword_symbols[j];
                    codeword_amps[i] = codeword_amps[j];
                }

                /* newest symbols at end of buffer (uses final i from last loop) */
                for(i=(interleave_frames-1)*coded_syms_per_frame,j=0; 
                    i<interleave_frames*coded_syms_per_frame; 
                    i++,j++) {
                    codeword_symbols[i] = payload_syms[j];
                    codeword_amps[i]    = payload_amps[j];
                }

                /* run de-interleaver */
                COMP  codeword_symbols_de[interleave_frames*coded_syms_per_frame];
                float codeword_amps_de[interleave_frames*coded_syms_per_frame];

                gp_deinterleave_comp (codeword_symbols_de, codeword_symbols, 
                                        interleave_frames*coded_syms_per_frame);
                gp_deinterleave_float(codeword_amps_de, codeword_amps, 
                                        interleave_frames*coded_syms_per_frame);

                float llr[coded_bits_per_frame];

                if (config_ldpc_en) {
                    uint8_t out_char[coded_bits_per_frame];

                    interleaver_sync_state_machine(ofdm, &ldpc, ofdm_config, 
                                codeword_symbols_de, codeword_amps_de, EsNo,
                                interleave_frames, iter, parityCheckCount, Nerrs_coded);

                    if ((ofdm->sync_state_interleaver == synced) && 
                            (ofdm->frame_count_interleaver == interleave_frames)) {
                        ofdm->frame_count_interleaver = 0;

                        if (config_testframes) {
                            Terrs += count_uncoded_errors(&ldpc, ofdm_config, Nerrs_raw, 
                                        interleave_frames, codeword_symbols_de);
                            Tbits += coded_bits_per_frame*interleave_frames; 
                        }

                        for (j=0; j<interleave_frames; j++) {
                            symbols_to_llrs(llr, &codeword_symbols_de[j*coded_syms_per_frame],
                                                 &codeword_amps_de[j*coded_syms_per_frame],
                                                 EsNo, ofdm->mean_amp, coded_syms_per_frame);
                            iter[j] = run_ldpc_decoder(&ldpc, out_char, llr, &parityCheckCount[j]);

                            //fprintf(stderr,"j: %d iter: %d pcc: %d\n", j, iter[j], parityCheckCount[j]);

                            if (config_testframes) {
                                /* construct payload data bits */
                                uint8_t payload_data_bits[data_bits_per_frame];
                                ofdm_generate_payload_data_bits(payload_data_bits, data_bits_per_frame);
                                
                                Nerrs_coded[j] = count_errors(payload_data_bits, out_char, data_bits_per_frame);
                                Terrs_coded += Nerrs_coded[j];
                                Tbits_coded += data_bits_per_frame;
                            }

                            fwrite(out_char, sizeof(char), data_bits_per_frame, fout);
                        }
                    } /* if interleaver synced ..... */
                } else { 
                    /* lpdc_en == 0,  external LDPC decoder, so output LLRs */
                    symbols_to_llrs(llr, codeword_symbols_de, codeword_amps_de, EsNo, ofdm->mean_amp, coded_syms_per_frame);
                    fwrite(llr, sizeof(double), coded_bits_per_frame, fout);
                }
            } else {    // !llrs_en (or ldpc_en)

                /* simple hard decision output for uncoded testing, all bits in frame dumped inlcuding UW and txt */
                for(i=0; i<ofdm_bitsperframe; i++) {
                    rx_bits_char[i] = rx_bits[i];
                }
                fwrite(rx_bits_char, sizeof(char), ofdm_bitsperframe, fout);
            }

            /* optional error counting on uncoded data in non-LDPC testframe mode */

            if (config_testframes && (config_ldpc_en == 0)) {
                /* build up a test frame consisting of unique word, txt bits, and psuedo-random
                   uncoded payload bits.  The psuedo-random generator is the same as Octave so
                   it can interoperate with ofdm_tx.m/ofdm_rx.m */

                int Npayloadbits = ofdm_bitsperframe-(ofdm_nuwbits+ofdm_ntxtbits);
                uint16_t r[Npayloadbits];
                uint8_t  payload_bits[Npayloadbits];
                uint8_t  tx_bits[Npayloadbits];
                
                ofdm_rand(r, Npayloadbits);

                for(i=0; i<Npayloadbits; i++) {
                    payload_bits[i] = r[i] > 16384;
                    //fprintf(stderr,"%d %d ", r[j], tx_bits_char[i]);
                }

                uint8_t txt_bits[ofdm_ntxtbits];

                for(i=0; i<ofdm_ntxtbits; i++) {
                    txt_bits[i] = 0;
                }

                ofdm_assemble_modem_frame(ofdm, tx_bits, payload_bits, txt_bits);

                Nerrs = 0;
                for(i=0; i<ofdm_bitsperframe; i++) {
                    if (tx_bits[i] != rx_bits[i]) {
                        Nerrs++;
                    }
                }

                Terrs += Nerrs;
                Tbits += ofdm_bitsperframe;

                if (frame_count >= NDISCARD) {
                    Terrs2 += Nerrs;
                    Tbits2 += ofdm_bitsperframe;
                }
            } // config_testframes ...

            frame_count++;
        } // state "synced" or "trial"

        nin_frame = ofdm_get_nin(ofdm);
        ofdm_sync_state_machine(ofdm, rx_uw);

        /* act on any events returned by state machine */

        if (ofdm->sync_start) {
            Terrs = Tbits = Terrs2 = Tbits2 = Terrs_coded = Tbits_coded = frame_count = 0;
            for(i=0; i<interleave_frames; i++) {
                Nerrs_raw[i] = Nerrs_coded[i] = 0;
            }

        }

        int r = 0;
        if (config_testframes && config_verbose) {
            r = (ofdm->frame_count_interleaver - 1 ) % interleave_frames;
            fprintf(stderr, "%3d st: %-6s", f, statemode[ofdm->last_sync_state]);
            fprintf(stderr, " euw: %2d %1d f: %5.1f ist: %-6s %2d eraw: %3d ecdd: %3d iter: %3d pcc: %3d",
                ofdm->uw_errors, ofdm->sync_counter,
                (double)ofdm->foff_est_hz,
                statemode[ofdm->last_sync_state_interleaver], 
                ofdm->frame_count_interleaver,
                Nerrs_raw[r], Nerrs_coded[r], iter[r], parityCheckCount[r]);
            fprintf(stderr, "\n");
        }

        if (config_log_payload_syms) {
            if (! log_payload_syms_flag) {
                memset(payload_syms, 0, (sizeof(COMP)*coded_syms_per_frame));
                memset(payload_amps, 0, (sizeof(float)*coded_syms_per_frame));
                }
            fwrite(payload_syms, sizeof(COMP), coded_syms_per_frame, fdiag);
            fwrite(payload_amps, sizeof(float), coded_syms_per_frame, fdiag);
            }

        f++;
    } // while(fread(.., fin))

    flush_all();    // To make sure this function is included in binary.
    fclose(fin);
    fclose(fout);
    fclose(fdiag);

    if (config_testframes) {
        printf("BER......: %5.4f Tbits: %5d Terrs: %5d\n", (double)Terrs/Tbits, Tbits, Terrs);
        if (!config_ldpc_en) {
            printf("BER2.....: %5.4f Tbits: %5d Terrs: %5d\n", (double)Terrs2/Tbits2, Tbits2, Terrs2);
        }
        if (config_ldpc_en) {
            printf("Coded BER: %5.4f Tbits: %5d Terrs: %5d\n",
                    (double)Terrs_coded/Tbits_coded, Tbits_coded, Terrs_coded);
        }
    }

    if (config_profile) {
        printf("\nStart Profile Data\n");
        machdep_profile_print_logged_samples();
        printf("End Profile Data\n");
        }

    printf("\nEnd of Test\n");
    fclose(stdout);
    fclose(stderr);

    return 0;
}

/* vi:set ts=4 et sts=4: */
