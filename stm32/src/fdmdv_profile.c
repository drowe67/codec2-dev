/*---------------------------------------------------------------------------*\

  FILE........: fdmdv_profile.c
  AUTHOR......: David Rowe
  DATE CREATED: 18 July 2014

  Profiling FDMDV modem operation on the STM32F4.

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
#include <stdint.h>
#include <math.h>

#include "stm32f4xx_conf.h"
#include "stm32f4xx.h"
#include "gdb_stdio.h"
#include "codec2_fdmdv.h"
#include "dump.h"
#include "sine.h"
#include "machdep.h"

#ifdef __EMBEDDED__
#define printf gdb_stdio_printf
#define fopen gdb_stdio_fopen
#define fclose gdb_stdio_fclose
#define fread gdb_stdio_fread
#define fwrite gdb_stdio_fwrite
#endif

#define TEST_FRAMES 25
#define CHANNEL_BUF_SIZE (10*FDMDV_NOM_SAMPLES_PER_FRAME)

static int  channel_count = 0;
static COMP channel[CHANNEL_BUF_SIZE];

static void channel_in(COMP tx_fdm[], int nout) {
    int i;

    /* add M_PITCH tx samples to end of buffer */

    assert((channel_count + nout) < CHANNEL_BUF_SIZE);
    for(i=0; i<nout; i++)
        channel[channel_count+i] = tx_fdm[i];
    channel_count += M_PITCH;
}

static void channel_out(COMP rx_fdm[], int nin) {
    int i,j;

    /* take nin samples from start of buffer */

    for(i=0; i<nin; i++) {
        rx_fdm[i] = channel[i];
    }

    /* shift buffer back */

    for(i=0,j=nin; j<channel_count; i++,j++)
        channel[i] = channel[j];
    channel_count -= nin;
}

int main(int argc, char *argv[]) {
    struct FDMDV       *fdmdv;
    int                 bits_per_fdmdv_frame,  bits_per_codec_frame;
    int                *tx_bits;
    int                *rx_bits;
    int                *codec_bits;
    COMP                tx_fdm[2*FDMDV_NOM_SAMPLES_PER_FRAME];
    COMP                rx_fdm[FDMDV_NOM_SAMPLES_PER_FRAME];
    int                 i, j, nin, reliable_sync_bit[2], sync_bit, bit_errors, ntest_bits, test_frame_sync;
    short              *error_pattern;
    struct MODEM_STATS  stats;
    PROFILE_VAR(mod_start, demod_start);

    machdep_profile_init ();
    fdmdv = fdmdv_create(FDMDV_NC);

    bits_per_fdmdv_frame = fdmdv_bits_per_frame(fdmdv);
    bits_per_codec_frame = 2*fdmdv_bits_per_frame(fdmdv);
    tx_bits = (int*)malloc(sizeof(int)*bits_per_codec_frame); assert(tx_bits != NULL);
    rx_bits = (int*)malloc(sizeof(int)*bits_per_codec_frame); assert(rx_bits != NULL);
    codec_bits = (int*)malloc(sizeof(int)*bits_per_codec_frame); assert(rx_bits != NULL);
    error_pattern = (short*)malloc(fdmdv_error_pattern_size(fdmdv)*sizeof(int)); assert(error_pattern != NULL);

    nin = FDMDV_NOM_SAMPLES_PER_FRAME;
    test_frame_sync = 0;

    for(i=0; i<TEST_FRAMES; i++) {
	fdmdv_get_test_bits(fdmdv, tx_bits);
	fdmdv_get_test_bits(fdmdv, &tx_bits[bits_per_fdmdv_frame]);

        PROFILE_SAMPLE(mod_start);

	fdmdv_mod(fdmdv, tx_fdm, tx_bits, &sync_bit);
	assert(sync_bit == 1);
	fdmdv_mod(fdmdv, &tx_fdm[FDMDV_NOM_SAMPLES_PER_FRAME], &tx_bits[bits_per_fdmdv_frame], &sync_bit);
	assert(sync_bit == 0);
        channel_in(tx_fdm, 2*FDMDV_NOM_SAMPLES_PER_FRAME);

        PROFILE_SAMPLE_AND_LOG(demod_start, mod_start, "  mod");

        for(j=0; j<2; j++) {
            channel_out(rx_fdm, nin);
            fdmdv_demod(fdmdv, rx_bits, &reliable_sync_bit[j], rx_fdm, &nin);
            if (reliable_sync_bit[j] == 0)
                memcpy(codec_bits, rx_bits, bits_per_fdmdv_frame*sizeof(int));
            else {
                memcpy(&codec_bits[bits_per_fdmdv_frame], rx_bits, bits_per_fdmdv_frame*sizeof(int));
                fdmdv_put_test_bits(fdmdv, &test_frame_sync, error_pattern, &bit_errors, &ntest_bits, codec_bits);
                fdmdv_put_test_bits(fdmdv, &test_frame_sync, error_pattern, &bit_errors, &ntest_bits, &codec_bits[bits_per_fdmdv_frame]);
            }
        }
        PROFILE_SAMPLE_AND_LOG2(demod_start, "  demod");
        PROFILE_SAMPLE_AND_LOG2(mod_start, "  mod & demod");

        fdmdv_get_demod_stats(fdmdv, &stats);

        printf("frame: %d sync: %d reliable_sync_bit: %d %d SNR: %3.2f test_frame_sync: %d\n",
               i, stats.sync, reliable_sync_bit[0], reliable_sync_bit[1], (double)stats.snr_est,
               test_frame_sync);
        machdep_profile_print_logged_samples();
    }

    fdmdv_destroy(fdmdv);

    return 0;
}

