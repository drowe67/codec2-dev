/*---------------------------------------------------------------------------*\

  FILE........: adc_sfdr_ut.c
  AUTHOR......: David Rowe
  DATE CREATED: August 2015

  Unit test for high speed ADC SFDR testing.  Samples ADC1 from in PA1 at
  Fs=2 MHz and write raw samples to a file, in discontinuus blocks of
  ADC_TUNER_BUF_SZ/2 samples.  The blocks are discontinuous as we
  don'thave the bandwitdh back to the host to support continuous sampling.

  To process the blocks, fread() ADC_TUNER_BUF_SZ/2 samples at a time,
  abs(fft) and sum results from next block.

\*---------------------------------------------------------------------------*/

/*
  Copyright (C) 2015 David Rowe

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
#include <stdlib.h>
#include "gdb_stdio.h"
#include "stm32f4_dac.h"
#include "stm32f4_adc_tuner.h"
#include "iir_tuner.h"
#include "sm1000_leds_switches.h"
#include "../src/codec2_fm.h"
#include "stm32f4xx.h"

#define BUFS          10
#define FS            2E6
#define N             1024

extern int adc_overflow1;

int main(void) {
    unsigned short unsigned_buf[N];
    short          buf[N];
    int            sam;
    int            i, j, fifo_sz;
    FILE          *fadc;

    fadc = fopen("adc.raw", "wb");
    if (fadc == NULL) {
        printf("Error opening output file: adc.raw\n\nTerminating....\n");
        exit(1);
    }
    fifo_sz = ADC_TUNER_BUF_SZ;
    printf("Starting! bufs: %d %d\n", BUFS, fifo_sz);

    adc_open(fifo_sz);
    adc_set_tuner_en(0); /* dump raw samples, no tuner */

    sm1000_leds_switches_init();

    for (i=0; i<BUFS; i++) {
        while(adc1_read((short*)unsigned_buf, N) == -1);

        /* convert to signed */

        for(j=0; j<N; j++) {
            sam = (int)unsigned_buf[j] - 32768;
            buf[j] = sam;
        }

        /* most of the time will be spent here */

        GPIOE->ODR |= (1 << 3);
        fwrite(buf, sizeof(short), N, fadc);
        GPIOE->ODR &= ~(1 << 3);
    }
    fclose(fadc);

    printf("Finished!\n");
}
