/*---------------------------------------------------------------------------*\

  FILE........: tst_ofdm_mod.c
  AUTHOR......: David Rowe, Don Reid
  DATE CREATED: 25 June 2018

  Test and profile OFDM modulation on the STM32F4.

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

/* This is a unit test implementation of the OFDM Mod function.
 *
 * Typical run:

    ofdm_gen_test_bits stm_in.raw 6 --rand 

    ofdm_mod stm_in.raw ref_mod_out.raw

    echo "00000000" > stm_cfg.txt

    <Load stm32 and run>

    compare_ints -s -b2 ref_mod_out.raw mod.raw

    ofdm_demod ref_mod_out.raw ref_ofdm_demod.raw --testframes
    ofdm_demod mod.raw stm_demod.raw --testframes

 * For LDPC use:

    ofdm_gen_test_bits stm_in.raw 6 --rand --ldpc

    ofdm_mod stm_in.raw ref_mod_out.raw --ldpc

    echo "00100000" > stm_cfg.txt

    <Load stm32 and run>

    ofdm_demod ref_mod_out.raw ref_ofdm_demod.raw --ldpc --testframes
    ofdm_demod mod.raw stm_demod.raw --ldpc --testframes

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
#include "interldpc.h"
#include "gp_interleaver.h"

#include "stm32f4xx_conf.h"
#include "stm32f4xx.h"
#include "machdep.h"

int main(int argc, char *argv[]) {
    struct OFDM *ofdm;
    FILE        *fcfg, *fin, *fout;
    struct LDPC  ldpc;

    // Test configuration, read from stm_cfg.txt
//    int          config_verbose;
//    int          config_testframes;
    int          config_ldpc_en;
//    int          config_log_payload_syms;
    int          config_profile;

    int          n_bpf, n_spf;
    int          frame = 0;
    int          i, j;

    static const int interleave_frames = 1; // No interleaving yet!

    semihosting_init();

    printf("OFDM_mod test and profile\n");

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
//    config_verbose = config[0] - '0';
//    config_testframes = config[1] - '0';
    config_ldpc_en = config[2] - '0';
//    config_log_payload_syms = config[3] - '0';
    config_profile = config[4] - '0';
    fclose(fcfg);

    PROFILE_VAR(ofdm_mod_start);
    if (config_profile) machdep_profile_init();

    struct OFDM_CONFIG *ofdm_config;
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

    ofdm_config->rs = (1.0f / ofdm_config->ts); /* Symbol Rate */

    /* config options can change here, non yet */

    ofdm = ofdm_create(ofdm_config);
    assert(ofdm != NULL);

    free(ofdm_config);

    /* Get a copy of the actual modem config */
    ofdm_config = ofdm_get_config_param();

    set_up_hra_112_112(&ldpc, ofdm_config);

    if (config_ldpc_en) {
        n_bpf =  interleave_frames * ldpc.data_bits_per_frame;
    } else {
        n_bpf = ofdm_get_bits_per_frame(ofdm);
    }

    n_spf = ofdm_get_samples_per_frame();
//    int ofdm_nuwbits = (ofdm_config->ns - 1) * ofdm_config->bps -
//                                                    ofdm_config->txtbits;
    int ofdm_ntxtbits =  ofdm_config->txtbits;

    uint8_t tx_bits_char[n_bpf];
    int16_t tx_scaled[n_spf];
    uint8_t txt_bits_char[ofdm_ntxtbits*interleave_frames];

    for(i=0; i< ofdm_ntxtbits*interleave_frames; i++) {
        txt_bits_char[i] = 0;
    }
   
    fin = fopen("stm_in.raw", "rb");
    if (fin == NULL) {
        printf("Error opening input file\n");
        exit(1);
    }

    fout = fopen("mod.raw", "wb");
    if (fout == NULL) {
        printf("Error opening output file\n");
        exit(1);
    }

    while (fread(tx_bits_char, sizeof(char), n_bpf, fin) == n_bpf) {

        if (config_profile) PROFILE_SAMPLE(ofdm_mod_start);

            if (config_ldpc_en) {

                complex float tx_sams[interleave_frames * n_spf];
                ofdm_ldpc_interleave_tx(ofdm, &ldpc, tx_sams, tx_bits_char,
                                txt_bits_char, interleave_frames, ofdm_config);

                for (j=0; j<interleave_frames; j++) {
                    for(i=0; i<n_spf; i++) {
                        tx_scaled[i] = OFDM_AMP_SCALE * crealf(tx_sams[j * n_spf + i]);
                    }

                }

             } else { // !config_ldpc_en
                int tx_bits[n_bpf];

                for(i=0; i<n_bpf; i++) {
                    tx_bits[i] = tx_bits_char[i];
                }

                COMP tx_sams[n_spf];
                ofdm_mod(ofdm, tx_sams, tx_bits);

                for(i=0; i<n_spf; i++) {
                    tx_scaled[i] = OFDM_AMP_SCALE * tx_sams[i].real;
                }
            }

        if (config_profile) PROFILE_SAMPLE_AND_LOG2(ofdm_mod_start, "  ofdm_mod");

        fwrite(tx_scaled, sizeof(int16_t), n_spf, fout);

        frame ++;

    }  // while (fread(...

    fclose(fin);
    fclose(fout);

    printf("%d frames processed\n", frame);

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
