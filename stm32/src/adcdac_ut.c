/*---------------------------------------------------------------------------*\

  FILE........: adcdac_ut.c
  AUTHOR......: David Rowe
  DATE CREATED: May 31 201310 Aug 2014

  Echoes ADC2 input (mic) to DAC2 output (speaker) on SM1000.

\*---------------------------------------------------------------------------*/

/*
  Copyright (C) 2013 David Rowe

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
#include "stm32f4_dac.h"
#include "stm32f4_adc.h"
#include "sm1000_leds_switches.h"

#define SINE_SAMPLES   32


/* 32 sample sine wave which at Fs=16kHz will be 500Hz.  Note samples
   are 16 bit 2's complement, the DAC driver convertsto 12 bit
   unsigned. */

short aSine[] = {
     -16,    6384,   12528,   18192,   23200,   27232,   30256,   32128,
   32752,   32128,   30256,   27232,   23152,   18192,   12528,    6384,
     -16,   -6416,  -12560,  -18224,  -23184,  -27264,  -30288,  -32160,
  -32768,  -32160,  -30288,  -27264,  -23184,  -18224,  -12560,   -6416
};

int main(void) {
    short buf[SINE_SAMPLES];
    int   i;

    dac_open(ADC_FS_16KHZ,4*DAC_BUF_SZ);
    adc_open(ADC_FS_16KHZ, 4*ADC_BUF_SZ);
    sm1000_leds_switches_init();

    while (1) {

        /* keep DAC FIFOs topped up */

        while(adc2_read(buf, SINE_SAMPLES) == -1);

        if (!switch_select()) {
            for(i=0; i<SINE_SAMPLES; i++)
                buf[i] = aSine[i];
        }

        dac2_write(buf, SINE_SAMPLES);
    }

}

