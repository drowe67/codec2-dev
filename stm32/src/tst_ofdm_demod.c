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

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <errno.h>

#include "codec2_ofdm.h"
#include "ofdm_internal.h"
#include "mpdecode_core.h"
#include "interldpc.h"
#include "gp_interleaver.h"
#include "test_bits_ofdm.h"

#include "stm32f4xx_conf.h"
#include "stm32f4xx.h"
#include "gdb_stdio.h"
#include "machdep.h"

#ifdef __EMBEDDED__
#define printf gdb_stdio_printf
#define fopen gdb_stdio_fopen
#define fclose gdb_stdio_fclose
#define fread gdb_stdio_fread
#define fwrite gdb_stdio_fwrite
#endif

#define NDISCARD 20

extern const int payload_data_bits[];
extern const int test_bits_ofdm[];

int main(int argc, char *argv[]) {
    struct OFDM *ofdm;
    FILE        *fin, *fout, *fdiag;
    int          nin_frame;

    int          i, f, interleave_frames;
    int          Nerrs, Terrs, Tbits, Terrs2, Tbits2, testframes, frame_count;
    int          ldpc_en, Tbits_coded, Terrs_coded;

    // For now;
    interleave_frames = 1;
    ldpc_en = 0;
    testframes = 1;

    int Nerrs_raw[interleave_frames];
    int Nerrs_coded[interleave_frames];
    int iter[interleave_frames];
    int parityCheckCount[interleave_frames];

    for(i=0; i<interleave_frames; i++) {
        Nerrs_raw[i] = Nerrs_coded[i] = iter[i] = parityCheckCount[i] = 0;
    }

    printf("OFDM_demod test and profile\n");

    PROFILE_VAR(ofdm_demod_start, ofdm_demod_sync_search,
                ofdm_demod_demod, ofdm_demod_diss, ofdm_demod_snr);

    machdep_profile_init();

    ofdm = ofdm_create(NULL);
    ofdm_set_verbose(ofdm, 1);

    int Nbitsperframe = ofdm_get_bits_per_frame(ofdm);
    //int Nsamperframe = ofdm_get_samples_per_frame();
    int Nmaxsamperframe = ofdm_get_max_samples_per_frame();
    int coded_syms_per_frame = ((OFDM_BITSPERFRAME/OFDM_BPS) -
                                (OFDM_NUWBITS/OFDM_BPS) -
                                (OFDM_NTXTBITS/OFDM_BPS));


    short  rx_scaled[Nmaxsamperframe];
    COMP   rxbuf_in[Nmaxsamperframe];
    int    rx_bits[Nbitsperframe];
    char   rx_bits_char[Nbitsperframe];
    int    rx_uw[OFDM_NUWBITS];
    short  txt_bits[OFDM_NTXTBITS];
    f = 0;
    Nerrs = Terrs = Tbits = Terrs2 = Tbits2 = Terrs_coded = Tbits_coded = frame_count = 0;

    float EsNo = 3;
    printf("Warning EsNo: %f hard coded\n", (double)EsNo);

    float snr_est_smoothed_dB = 0.0;

    COMP  payload_syms[coded_syms_per_frame];
    float payload_amps[coded_syms_per_frame];

printf("coded_syms_per_frame = %d\n", coded_syms_per_frame);
printf("sizeof(COMP) = %d\n", sizeof(COMP));
printf("payload_syms = %d\n", sizeof(COMP)*coded_syms_per_frame);
printf("sizeof(float) = %d\n", sizeof(float));
printf("payload_amps = %d\n", sizeof(float)*coded_syms_per_frame);

    fin = fopen("stm_in.raw", "rb");
    if (fin == NULL) {
        printf("Error opening input file\n");
        exit(1);
    }

    fout = fopen("stm_out.raw", "wb");
    if (fout == NULL) {
        printf("Error opening output file\n");
        exit(1);
    }

    fdiag = fopen("stm_diag.raw", "wb");
    if (fdiag == NULL) {
        printf("Error opening diag file\n");
        exit(1);
    }

    nin_frame = ofdm_get_nin(ofdm);
    while(fread(rx_scaled, sizeof(short), nin_frame, fin) == nin_frame) {

        int log_payload_syms = 0;

        PROFILE_SAMPLE(ofdm_demod_start);

	    /* scale and demod */
	    for(i=0; i<nin_frame; i++) {
	        rxbuf_in[i].real = (float)rx_scaled[i]/(OFDM_AMP_SCALE/2);
                rxbuf_in[i].imag = 0.0;
            }

        PROFILE_SAMPLE_AND_LOG2(ofdm_demod_start, "  ofdm_demod_start");

        if (strcmp(ofdm->sync_state,"search") == 0) {
            PROFILE_SAMPLE(ofdm_demod_sync_search);
            ofdm_sync_search(ofdm, rxbuf_in);
            PROFILE_SAMPLE_AND_LOG2(ofdm_demod_sync_search, "  ofdm_demod_sync_search");
        }

        if ((strcmp(ofdm->sync_state,"synced") == 0) ||
            (strcmp(ofdm->sync_state,"trial") == 0) ) {
            PROFILE_SAMPLE(ofdm_demod_demod);
            ofdm_demod(ofdm, rx_bits, rxbuf_in);
            PROFILE_SAMPLE_AND_LOG2(ofdm_demod_demod, "  ofdm_demod_demod");
            PROFILE_SAMPLE(ofdm_demod_diss);
            ofdm_disassemble_modem_frame(ofdm, rx_uw, payload_syms, payload_amps, txt_bits);
            PROFILE_SAMPLE_AND_LOG2(ofdm_demod_diss, "  ofdm_demod_diss");
            log_payload_syms = 1;

            /* SNR estimation and smoothing */
            PROFILE_SAMPLE(ofdm_demod_snr);

            float snr_est_dB = 10*log10((ofdm->sig_var/ofdm->noise_var)*OFDM_NC*OFDM_RS/3000);
            snr_est_smoothed_dB = 0.9*snr_est_smoothed_dB + 0.1*snr_est_dB;
            PROFILE_SAMPLE_AND_LOG2(ofdm_demod_snr, "  ofdm_demod_snr");

            /* simple hard decision output for uncoded testing, all bits in frame dumped inlcuding UW and txt */
            for(i=0; i<Nbitsperframe; i++) {
                rx_bits_char[i] = rx_bits[i];
            }
            fwrite(rx_bits_char, sizeof(char), Nbitsperframe, fout);

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
            } // testframes ...

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

        int  r=0;
        if (testframes) {
            r = (ofdm->frame_count_interleaver - 1 ) % interleave_frames;
        }
        printf("%3d st: %-6s", f, ofdm->last_sync_state);
        printf(" euw: %2d %1d f: %5.1f ist: %-6s %2d eraw: %3d ecdd: %3d iter: %3d pcc: %3d",
                ofdm->uw_errors, ofdm->sync_counter,
                (double)ofdm->foff_est_hz,
                ofdm->last_sync_state_interleaver, ofdm->frame_count_interleaver,
                Nerrs_raw[r], Nerrs_coded[r], iter[r], parityCheckCount[r]);
        printf("\n");

        if (! log_payload_syms) {
            memset(payload_syms, 0, (sizeof(COMP)*coded_syms_per_frame));
            memset(payload_amps, 0, (sizeof(float)*coded_syms_per_frame));
            }
        fwrite(payload_syms, sizeof(COMP), coded_syms_per_frame, fdiag);
        fwrite(payload_amps, sizeof(float), coded_syms_per_frame, fdiag);

        f++;
    } // while(fread(.., fin))

    fclose(fin);
    fclose(fout);

    if (testframes) {
        printf("BER......: %5.4f Tbits: %5d Terrs: %5d\n", (double)Terrs/Tbits, Tbits, Terrs);
        if (!ldpc_en) {
            printf("BER2.....: %5.4f Tbits: %5d Terrs: %5d\n", (double)Terrs2/Tbits2, Tbits2, Terrs2);
        }
        if (ldpc_en) {
            printf("Coded BER: %5.4f Tbits: %5d Terrs: %5d\n",
                    (double)Terrs_coded/Tbits_coded, Tbits_coded, Terrs_coded);
        }        
    }

    printf("End of test\n");
    
    machdep_profile_print_logged_samples();

    return 0;
}

/* vi:set ts=4 et sts=4: */
