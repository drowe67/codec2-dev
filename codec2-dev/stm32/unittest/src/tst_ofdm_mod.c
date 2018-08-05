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

    ofdm_get_test_bits stm_in.raw -f 10

    ofdm_mod stm_in.raw ref_mod_out.raw

    <Load stm32 and run>

    compare_ints -s 2 ref_mod_out.raw mod.raw

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
    FILE        *fin, *fout;
    int          n_bpf, n_spf;
    int          i;

    semihosting_init();

    printf("OFDM_mod test and profile\n");

    PROFILE_VAR(ofdm_mod_start);

    machdep_profile_init();

    struct OFDM_CONFIG *ofdm_config;
    if ((ofdm_config = (struct OFDM_CONFIG *) calloc(1, sizeof (struct OFDM_CONFIG))) == NULL) {
        fprintf(stderr, "Out of Memory\n");
        exit(1);
    }
    ofdm = ofdm_create(ofdm_config);
    assert(ofdm != NULL);

    free(ofdm_config);

    n_bpf = ofdm_get_bits_per_frame(ofdm);
    n_spf = ofdm_get_samples_per_frame();
    uint8_t tx_bits_char[n_bpf];
    short   tx_scaled[n_spf];

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
        PROFILE_SAMPLE(ofdm_mod_start);

        int tx_bits[n_bpf];
           for(i=0; i<n_bpf; i++)
                tx_bits[i] = tx_bits_char[i];
            COMP tx_sams[n_spf];
            ofdm_mod(ofdm, tx_sams, tx_bits);

        for(i=0; i<n_spf; i++)
                tx_scaled[i] = OFDM_AMP_SCALE * tx_sams[i].real;

        PROFILE_SAMPLE_AND_LOG2(ofdm_mod_start, "  ofdm_mod");

        fwrite(tx_scaled, sizeof(short), n_spf, fout);
    }
    fclose(fin);
    fclose(fout);
    
    fflush(stdout);
    stdout = freopen("stm_profile", "w", stdout);
    machdep_profile_print_logged_samples();

    fclose(stdout);
    fclose(stderr);

    return 0;
}

/* vi:set ts=4 et sts=4: */
