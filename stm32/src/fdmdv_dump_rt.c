/*---------------------------------------------------------------------------*\

  FILE........: fdmdv_dump_rt.c
  AUTHOR......: David Rowe
  DATE CREATED: 9 Sep 2014

 Runs the fdmdv demod in real time for a few seconds then dumps some
 modem info to a text file for plotting in Octave.  Way to verify the
 "from radio" SM1000 hardware, ADC, and demod on the SM1000.

 Requires FreeDV signal to be sent to CN6 of SM1000.

 Octave:

   load scatter.txt
   l=length(scatter)
   plot(scatter(:,1:2:l),scatter(:,2:2:l),'+')

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

#include <stm32f4xx_gpio.h>
#include "stm32f4_adc.h"
#include "stm32f4_dac.h"
#include "freedv_api.h"
#include "codec2_fdmdv.h"
#include "sm1000_leds_switches.h"
#include "gdb_stdio.h"

#ifdef __EMBEDDED__
#define printf gdb_stdio_printf
#define fprintf gdb_stdio_fprintf
#define fopen gdb_stdio_fopen
#define fclose gdb_stdio_fclose
#define fread gdb_stdio_fread
#define fwrite gdb_stdio_fwrite
#endif

#define FREEDV_NSAMPLES_16K (2*FREEDV_NSAMPLES)
#define START_LOG_FRAMES 100
#define LOG_FRAMES       10
#define STOP_LOG_FRAMES  (START_LOG_FRAMES+LOG_FRAMES)
#define NC  16

int main(void) {
    struct freedv *f;
    short          adc16k[FDMDV_OS_TAPS_16K+FREEDV_NSAMPLES_16K];
    short          dac16k[FREEDV_NSAMPLES_16K];
    short          adc8k[FREEDV_NSAMPLES];
    short          dac8k[FDMDV_OS_TAPS_8K+FREEDV_NSAMPLES];

    int            nin, nout, i, j, frames, lines;

    COMP          *symb, *psymb;

    /* init all the drivers for various peripherals */

    sm1000_leds_switches_init();
    dac_open(4*DAC_BUF_SZ);
    adc_open(4*ADC_BUF_SZ);
    f = freedv_open(FREEDV_MODE_1600);

    /* clear filter memories */

    for(i=0; i<FDMDV_OS_TAPS_16K; i++)
	adc16k[i] = 0.0;
    for(i=0; i<FDMDV_OS_TAPS_8K; i++)
	dac8k[i] = 0.0;

    /* allocate storage for the symbols */

#define TMP
#ifdef TMP
    symb = (COMP*)malloc(sizeof(COMP)*(NC+1)*(STOP_LOG_FRAMES - START_LOG_FRAMES));
    assert(symb != NULL);
    psymb = symb;
    frames = 0;
    lines = 0;
#endif
    while(1)  {

        /* Receive --------------------------------------------------------------------------*/

        /* ADC1 is the demod in signal from the radio rx, DAC2 is the SM1000 speaker */

        nin = freedv_nin(f);
        nout = nin;
        f->total_bit_errors = 0;

        if (adc1_read(&adc16k[FDMDV_OS_TAPS_16K], 2*nin) == 0) {
            GPIOE->ODR = (1 << 3);
            fdmdv_16_to_8_short(adc8k, &adc16k[FDMDV_OS_TAPS_16K], nin);
            nout = freedv_rx(f, &dac8k[FDMDV_OS_TAPS_8K], adc8k);
            fdmdv_8_to_16_short(dac16k, &dac8k[FDMDV_OS_TAPS_8K], nout);
            dac2_write(dac16k, 2*nout);
            led_ptt(0); led_rt(f->fdmdv_stats.sync); led_err(f->total_bit_errors);
            GPIOE->ODR &= ~(1 << 3);

#define TMP1
#ifdef TMP1
            if (f->fdmdv_stats.sync)
                frames++;
            if ((frames >= START_LOG_FRAMES) && (lines < LOG_FRAMES)) {
                for(i=0; i<=f->fdmdv_stats.Nc; i++)
                    psymb[i] = f->fdmdv_stats.rx_symbols[i];
                psymb += (f->fdmdv_stats.Nc+1);
                lines++;
            }

            if (frames >= STOP_LOG_FRAMES) {
                FILE *ft = fopen("scatter.txt", "wt");
                assert(ft != NULL);
                printf("Writing scatter file....\n");
                for(j=0; j<LOG_FRAMES; j++) {
                    for(i=0; i<=f->fdmdv_stats.Nc; i++) {
                        fprintf(ft, "%f\t%f\t",
                                (double)symb[j*(f->fdmdv_stats.Nc+1)+i].real,
                                (double)symb[j*(f->fdmdv_stats.Nc+1)+i].imag);
                        printf("line: %d\n", j);
                    }
                    fprintf(ft, "\n");
                }
                fclose(ft);
                printf("SNR = %3.2f dB\nfinished!\n", (double)f->fdmdv_stats.snr_est);
                while(1);
            }
#endif
        }

    } /* while(1) ... */
}

