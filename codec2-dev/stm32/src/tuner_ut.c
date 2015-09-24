/*---------------------------------------------------------------------------*\

  FILE........: tuner_ut.c
  AUTHOR......: David Rowe
  DATE CREATED: 20 Feb 2015

  Unit test for high speed ADC radio tuner, samples signal centred at
  500kHz using Fs=2 MHz and uploads to host at Fs=10 kHz.

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

#define REC_TIME_SECS 10
#define FS            50000
#define N             5000

extern int adc_overflow1;

int main(void) {
    float  tuner_out[IIR_TUNER_DEC_50_10_FILT_MEM+N/2];
    float  fm_out[N/2];
    //float  dec_10[(N/2)/5];
    short  dec_10_short[(N/2)];
    int    bufs, i, j, k, fifo_sz, bn;
    FILE  *ftuner;
    struct FM *fm;

    ftuner = fopen("tuner.raw", "wb");
    if (ftuner == NULL) {
        printf("Error opening input file: tuner.raw\n\nTerminating....\n");
        exit(1);
    }
    bufs = FS*REC_TIME_SECS/N;
    fifo_sz = ((4*N/ADC_TUNER_N)+1)*ADC_TUNER_N;
    printf("Starting! bufs: %d %d\n", bufs, fifo_sz);

    //dac_open(DAC_BUF_SZ);
    adc_open(fifo_sz);
    sm1000_leds_switches_init();

    fm         = fm_create(N/2);
    fm->Fs     = 44400.0;
    fm->fm_max = 3000.0;
    fm->fd     = 5000.0;
    fm->fc     = fm->Fs/4;

    i = 0; bn = 0;
    while(1) {
        /* wait for buffer of Fs=50kHz tuner output samples */

        while(adc1_read((short *)&tuner_out[IIR_TUNER_DEC_50_10_FILT_MEM], N) == -1);

        /* The semi-hosting system can only handle Fs=16kHz and below so resample down
           to Fs=10 kHz and convert to shorts */

        #ifdef SSB
        iir_tuner_dec_50_to_10(dec_10, &tuner_out[IIR_TUNER_DEC_50_10_FILT_MEM], N/2);
        for(j=0; j<IIR_TUNER_DEC_50_10_FILT_MEM; j++)
            tuner_out[j] = tuner_out[j+N/2];
        for(j=0; j<(N/2)/5; j++)
            dec_10_short[j] = dec_10[j]/ADC_TUNER_M;
        #else
        GPIOE->ODR |= (1 << 3);
        /*
        for(j=0; j<N/2; j++)
            tuner_out[IIR_TUNER_DEC_50_10_FILT_MEM+j] = 0;
        for(j=1; j<N/2; j+=4)
            tuner_out[IIR_TUNER_DEC_50_10_FILT_MEM+j] = 100;
        for(j=3; j<N/2; j+=4)
            tuner_out[IIR_TUNER_DEC_50_10_FILT_MEM+j] = -100;
        */
        fm_demod(fm, fm_out, &tuner_out[IIR_TUNER_DEC_50_10_FILT_MEM]);
        GPIOE->ODR &= ~(1 << 3);

        for(j=0,k=0; j<N/2; j+=5,k++)
            dec_10_short[k+bn] = 16384*fm_out[j];
        bn += (N/2)/5;
        #endif

        if ((bn == N/2) && (i < bufs)) {
            //for(j=0; j<N/2; j++)
            //    dec_10_short[j] = tuner_out[IIR_TUNER_DEC_50_10_FILT_MEM+j];
            fwrite(dec_10_short, sizeof(short), (N/2), ftuner);
            bn = 0;
            i += 5;
            //fwrite(dec_10_short, sizeof(short), (N/2)/5, ftuner);
            //printf("%d %d\n", i, adc_overflow1);
        }
        if (i == bufs) {
            printf("finished! %d\n", adc_overflow1);
            fclose(ftuner);
            i++;
        }
    }
}

