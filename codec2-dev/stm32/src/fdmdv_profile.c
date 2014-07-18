/*---------------------------------------------------------------------------*\

  FILE........: fdmdv_profile.c
  AUTHOR......: David Rowe
  DATE CREATED: 18 July 2014

  Profiling Codec 2 operation on the STM32F4.

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

int main(int argc, char *argv[]) {
    struct FDMDV       *fdmdv;
    int                 bits_per_fdmdv_frame,  bits_per_codec_frame;
    int                *tx_bits;
    int                *rx_bits;
    COMP                tx_fdm[2*FDMDV_NOM_SAMPLES_PER_FRAME];
    int                 i, nin, reliable_sync_bit, sync_bit;
    struct FDMDV_STATS  stats;
    TIMER_VAR(mod_start, demod_start);

    fdmdv = fdmdv_create(FDMDV_NC);

    bits_per_fdmdv_frame = fdmdv_bits_per_frame(fdmdv);
    bits_per_codec_frame = 2*fdmdv_bits_per_frame(fdmdv);
    tx_bits = (int*)malloc(sizeof(int)*bits_per_codec_frame); assert(tx_bits != NULL);
    rx_bits = (int*)malloc(sizeof(int)*bits_per_codec_frame); assert(rx_bits != NULL);

    nin = FDMDV_NOM_SAMPLES_PER_FRAME;

    for(i=0; i<TEST_FRAMES; i++) {
	fdmdv_get_test_bits(fdmdv, tx_bits);
	fdmdv_get_test_bits(fdmdv, &tx_bits[bits_per_fdmdv_frame]);

        TIMER_SAMPLE(mod_start);

	fdmdv_mod(fdmdv, tx_fdm, tx_bits, &sync_bit);
	assert(sync_bit == 1);
	fdmdv_mod(fdmdv, &tx_fdm[FDMDV_NOM_SAMPLES_PER_FRAME], &tx_bits[bits_per_fdmdv_frame], &sync_bit);
	assert(sync_bit == 0);

        TIMER_SAMPLE_AND_LOG(demod_start, mod_start, "  enc");     

        fdmdv_demod(fdmdv, rx_bits, &reliable_sync_bit, tx_fdm, &nin);
        fdmdv_demod(fdmdv, rx_bits, &reliable_sync_bit, &tx_fdm[FDMDV_NOM_SAMPLES_PER_FRAME], &nin);
        TIMER_SAMPLE_AND_LOG2(demod_start, "  dec");     
        TIMER_SAMPLE_AND_LOG2(mod_start, "  enc & dec");     

        fdmdv_get_demod_stats(fdmdv, &stats);
        printf("frame: %d sync: %d reliable_sync_bit: %d SNR: %3.2f\n", i, stats.sync, reliable_sync_bit, (double)stats.snr_est);
        machdep_timer_print_logged_samples();
    }

    fdmdv_destroy(fdmdv);

    return 0;
}

