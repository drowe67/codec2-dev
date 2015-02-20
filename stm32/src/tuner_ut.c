/*---------------------------------------------------------------------------*\

  FILE........: tuner_ut.c
  AUTHOR......: David Rowe
  DATE CREATED: 20 Feb 2015

  Unit test for high speed ADC radio tuner.

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
#include "stm32f4_dac.h"
#include "stm32f4_adc_tuner.h"
#include "sm1000_leds_switches.h"

int main(void) {
    short buf[SINE_SAMPLES];
    int   i;

    dac_open(4*DAC_BUF_SZ);
    adc_open(4*ADC_TUNER_N;
    sm1000_leds_switches_init();

    while (1) {

        while(adc1_read(buf, ADC_TUNER_N) == -1);
        
        /*
        if (switch_select()) {
            for(i=0; i<SINE_SAMPLES; i++)
                buf[i] = 0;
        }
            
        dac1_write(buf, SINE_SAMPLES);
        */
    }
   
}

