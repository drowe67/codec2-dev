/*---------------------------------------------------------------------------*\

  FILE........: freedv_rx_profile.c
  AUTHOR......: David Rowe
  DATE CREATED: 13 August 2014

  Profiling freedv_rx() operation on the STM32F4.

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

#define PROFILE

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>

#include "stm32f4xx_conf.h"
#include "stm32f4xx.h"
#include "gdb_stdio.h"
#include "freedv_api.h"
#include "machdep.h"
#include "codec2_fdmdv.h"

#ifdef __EMBEDDED__
#define printf gdb_stdio_printf
#define fopen gdb_stdio_fopen
#define fclose gdb_stdio_fclose
#define fread gdb_stdio_fread
#define fwrite gdb_stdio_fwrite
#define fprintf gdb_stdio_fprintf
#endif

#define FREEDV_NSAMPLES_16K (2*FREEDV_NSAMPLES)

int main(int argc, char *argv[]) {
    struct freedv *f;
    FILE          *fin, *fout, *ftotal;
    int            frame, nin_16k, nin, i, nout = 0;
    int            n_samples, n_samples_16k;
    int            sync;
    float          snr_est;

    PROFILE_VAR(fdmdv_16_to_8_start, freedv_rx_start, fdmdv_8_to_16_start);

    machdep_profile_init();

    f = freedv_open(FREEDV_MODE_1600);
    n_samples = freedv_get_n_speech_samples(f);
    n_samples_16k = 2*n_samples;

    short          adc16k[FDMDV_OS_TAPS_16K+n_samples_16k];
    short          dac16k[n_samples_16k];
    short          adc8k[n_samples];
    short          dac8k[FDMDV_OS_TAPS_8K+n_samples];

    // Receive ---------------------------------------------------------------------

    frame = 0;

    fin = fopen("mod_16k.raw", "rb");
    if (fin == NULL) {
        printf("Error opening input file\n");
        exit(1);
    }

    fout = fopen("speechout_16k.raw", "wb");
    if (fout == NULL) {
        printf("Error opening output file\n");
        exit(1);
    }

    ftotal = fopen("total.txt", "wt");
    assert(ftotal != NULL);

    /* clear filter memories */

    for(i=0; i<FDMDV_OS_TAPS_16K; i++)
	adc16k[i] = 0.0;
    for(i=0; i<FDMDV_OS_TAPS_8K; i++)
	dac8k[i] = 0.0;

    nin = freedv_nin(f);
    nin_16k = 2*nin;
    nout = nin;
    while (fread(&adc16k[FDMDV_OS_TAPS_16K], sizeof(short), nin_16k, fin) == nin_16k) {

        PROFILE_SAMPLE(fdmdv_16_to_8_start);

        fdmdv_16_to_8_short(adc8k, &adc16k[FDMDV_OS_TAPS_16K], nin);

        PROFILE_SAMPLE_AND_LOG(freedv_rx_start, fdmdv_16_to_8_start, "  fdmdv_16_to_8");

        nout = freedv_rx(f, &dac8k[FDMDV_OS_TAPS_8K], adc8k);
        nin = freedv_nin(f); nin_16k = 2*nin;

        PROFILE_SAMPLE_AND_LOG(fdmdv_8_to_16_start, freedv_rx_start, "  freedv_rx");

        fdmdv_8_to_16_short(dac16k, &dac8k[FDMDV_OS_TAPS_8K], nout);

        PROFILE_SAMPLE_AND_LOG2(fdmdv_8_to_16_start, "  fdmdv_8_to_16");

        fprintf(ftotal, "%d\n", machdep_profile_sample() - fdmdv_16_to_8_start);
        machdep_profile_print_logged_samples();

        fwrite(dac16k, sizeof(short), 2*nout, fout);
        freedv_get_modem_stats(f, &sync, &snr_est);
        printf("frame: %d nin_16k: %d sync: %d SNR: %3.2f \n",
               ++frame, nin_16k, sync, (double)snr_est);
    }

    fclose(fin);
    fclose(fout);
    fclose(ftotal);

    return 0;
}

