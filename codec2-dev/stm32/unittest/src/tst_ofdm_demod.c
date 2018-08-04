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
 *  tst_ofdm_demod_ideal
 *  tst_ofdm_demod_AWGN
 *  tst_ofdm_demod_fade
 *
 * See tst_ofdm_demod_setup and tst_ofdm_demod_check scripts for details.
 *
 * This program reads a file "stm_cfg.txt" at startup to configure its options.
 *
 * This program is intended to be run using input data, typically
 * Codec2 frames, which may have had simulated RF degredation applied.
 * For example:
 *
 *    ofdm_get_test_bits - 10 | * ofdm_mod - - | cohpsk_ch - stm_in.raw -20 -Fs 8000 -f -5
 *
 * Reference data can be created by running the same input through the x86 ofdm_demod tool.
 *
 *    ofdm_demod stm_in.raw ref_demod_out.raw -o ofdm_demod_ref_log.txt --testframes
 *
 * Comparison of the results to the reference will depend on the test conditions.
 *
 */


#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <errno.h>

#include "semihosting.h"
#include "codec2_ofdm.h"
#include "ofdm_internal.h"
#include "mpdecode_core.h"
#include "interldpc.h"
#include "gp_interleaver.h"
#include "test_bits_ofdm.h"

#include "stm32f4xx_conf.h"
#include "stm32f4xx.h"
#include "machdep.h"

#define NDISCARD 20

extern const int payload_data_bits[];
extern const int test_bits_ofdm[];

static struct OFDM_CONFIG *ofdm_config;

static int ofdm_bitsperframe;
static int ofdm_rowsperframe;
static int ofdm_nuwbits;
static int ofdm_ntxtbits;
static int ofdm_nin;

int main(int argc, char *argv[]) {
    struct OFDM *ofdm;
    FILE        *fcfg, *fin, *fout, *fdiag;
    int          nin_frame;

    int          config_verbose, config_testframes, config_ldpc_en, config_log_payload_syms,
                 config_profile;

    int          i, f;
    int          interleave_frames;
    int          Nerrs, Terrs, Tbits, Terrs2, Tbits2, frame_count;
    int          Tbits_coded, Terrs_coded;

    semihosting_init();

    fprintf(stdout, "OFDM Demod test\n");
    fprintf(stderr, "OFDM Demod test\n");

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
    if (config_profile) machdep_profile_init();

    if ((ofdm_config = (struct OFDM_CONFIG *) calloc(1, sizeof (struct OFDM_CONFIG))) == NULL) {
        fprintf(stderr, "Out of Memory\n");
        exit(1);
    }

    ofdm = ofdm_create(ofdm_config);
    assert(ofdm != NULL);

    free(ofdm_config);

    /* Get a copy of the actual modem config */
    ofdm_config = ofdm_get_config_param();

    ofdm_bitsperframe = ofdm_get_bits_per_frame();
    ofdm_rowsperframe = ofdm_bitsperframe / (ofdm_config->nc * ofdm_config->bps);
    ofdm_nuwbits = (ofdm_config->ns - 1) * ofdm_config->bps - ofdm_config->txtbits;
    ofdm_ntxtbits = ofdm_config->txtbits;
    ofdm_nin = ofdm_get_nin(ofdm);

    ofdm_set_verbose(ofdm, config_verbose);

    int Nmaxsamperframe = ofdm_get_max_samples_per_frame();
    int coded_syms_per_frame = ((ofdm_bitsperframe/ofdm_config->bps) -
                               (ofdm_nuwbits/ofdm_config->bps) -
                               (ofdm_ntxtbits/ofdm_config->bps));

    short  rx_scaled[Nmaxsamperframe];
    COMP   rxbuf_in[Nmaxsamperframe];
    int    rx_bits[ofdm_bitsperframe];
    char   rx_bits_char[ofdm_bitsperframe];
    int    rx_uw[ofdm_nuwbits];
    short  txt_bits[ofdm_ntxtbits];
    f = 0;
    Nerrs = Terrs = Tbits = Terrs2 = Tbits2 = Terrs_coded = Tbits_coded = frame_count = 0;

    float snr_est_smoothed_dB = 0.0;

    COMP  payload_syms[coded_syms_per_frame];
    float payload_amps[coded_syms_per_frame];

    fin = fopen("stm_in.raw", "rb");
    if (fin == NULL) {
        fprintf(stderr, "Error opening input file\n");
        exit(1);
    }

    fout = fopen("stm_out.raw", "wb");
    if (fout == NULL) {
        fprintf(stderr, "Error opening output file\n");
        exit(1);
    }

    fdiag = fopen("stm_diag.raw", "wb");
    if (fdiag == NULL) {
        fprintf(stderr, "Error opening diag file\n");
        exit(1);
    }

    nin_frame = ofdm_get_nin(ofdm);
    int num_read;
    while((num_read = fread(rx_scaled, sizeof(short), nin_frame, fin)) == nin_frame) {

        int log_payload_syms_flag = 0;

        if (config_profile) PROFILE_SAMPLE(ofdm_demod_start);

	    /* scale and demod */
	    for(i=0; i<nin_frame; i++) {
	        rxbuf_in[i].real = (float)rx_scaled[i]/(OFDM_AMP_SCALE/2);
                rxbuf_in[i].imag = 0.0;
            }

        if (config_profile) PROFILE_SAMPLE_AND_LOG2(ofdm_demod_start, "  ofdm_demod_start");

        if (strcmp(ofdm->sync_state,"search") == 0) {
            if (config_profile) PROFILE_SAMPLE(ofdm_demod_sync_search);
            ofdm_sync_search(ofdm, rxbuf_in);
            if (config_profile) PROFILE_SAMPLE_AND_LOG2(ofdm_demod_sync_search, "  ofdm_demod_sync_search");
        }

        if ((strcmp(ofdm->sync_state,"synced") == 0) ||
            (strcmp(ofdm->sync_state,"trial") == 0) ) {
            if (config_profile) PROFILE_SAMPLE(ofdm_demod_demod);
            ofdm_demod(ofdm, rx_bits, rxbuf_in);
            if (config_profile) PROFILE_SAMPLE_AND_LOG2(ofdm_demod_demod, "  ofdm_demod_demod");
            if (config_profile) PROFILE_SAMPLE(ofdm_demod_diss);
            ofdm_disassemble_modem_frame(ofdm, rx_uw, payload_syms, payload_amps, txt_bits);
            if (config_profile) PROFILE_SAMPLE_AND_LOG2(ofdm_demod_diss, "  ofdm_demod_diss");
            log_payload_syms_flag = 1;

            /* SNR estimation and smoothing */
            if (config_profile) PROFILE_SAMPLE(ofdm_demod_snr);

            float snr_est_dB = 10*log10((ofdm->sig_var/ofdm->noise_var) * ofdm_config->nc * ofdm_config->rs / 3000);
            snr_est_smoothed_dB = 0.9*snr_est_smoothed_dB + 0.1*snr_est_dB;
            if (config_profile) PROFILE_SAMPLE_AND_LOG2(ofdm_demod_snr, "  ofdm_demod_snr");

            /* simple hard decision output for uncoded testing, all bits in frame dumped inlcuding UW and txt */
            for(i=0; i<ofdm_bitsperframe; i++) {
                rx_bits_char[i] = rx_bits[i];
            }
            fwrite(rx_bits_char, sizeof(char), ofdm_bitsperframe, fout);

            /* optional error counting on uncoded data in non-LDPC testframe mode */

            if (config_testframes && (config_ldpc_en == 0)) {
                Nerrs = 0;
                for(i=0; i<ofdm_bitsperframe; i++) {
                    if (test_bits_ofdm[i] != rx_bits[i]) {
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
            fprintf(stderr, "%3d st: %-6s", f, ofdm->last_sync_state);
            fprintf(stderr, " euw: %2d %1d f: %5.1f ist: %-6s %2d eraw: %3d ecdd: %3d iter: %3d pcc: %3d",
                ofdm->uw_errors, ofdm->sync_counter,
                (double)ofdm->foff_est_hz,
                ofdm->last_sync_state_interleaver, ofdm->frame_count_interleaver,
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
        fflush(stdout);
        stdout = freopen("stm_profile", "w", stdout);
        machdep_profile_print_logged_samples();
        }

    fflush(stdout);
    fflush(stderr);

    return 0;
}

/* vi:set ts=4 et sts=4: */
