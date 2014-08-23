/*---------------------------------------------------------------------------*\

  FILE........: sm1000_main.c
  AUTHOR......: David Rowe
  DATE CREATED: August 5 2014

  Main program for SM1000.

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

#include "stm32f4_adc.h"
#include "stm32f4_dac.h"
#include "freedv_api.h"
#include "codec2_fdmdv.h"
#include "sm1000_leds_switches.h"
#include <stm32f4xx_gpio.h>

#define FREEDV_NSAMPLES_16K (2*FREEDV_NSAMPLES)

int main(void) {
    struct freedv *f;
    short          adc16k[FDMDV_OS_TAPS_16K+FREEDV_NSAMPLES_16K];
    short          dac16k[FREEDV_NSAMPLES_16K];
    short          adc8k[FREEDV_NSAMPLES];
    short          dac8k[FDMDV_OS_TAPS_8K+FREEDV_NSAMPLES];
 
    int    nin, nout, i;

    /* init all the drivers for various peripherals */

    //sm1000_leds_switches_init();
    dac_open(4*DAC_BUF_SZ);
    adc_open(4*ADC_BUF_SZ);
    f = freedv_open(FREEDV_MODE_1600);

    /* LEDs into a known state */

    //led_pwr(1); led_ptt(0); led_rt(0); led_err(0);

    /* 
       TODO:
       [ ] UT analog interfaces from file IO
       [ ] UTs for simultaneous tx & rx on analog interfaces
       [ ] measure CPU load of various parts with a blinky
           [ ] ADC and DAC drivers
           [ ] rate conversion
       [ ] detect program assert type errors with a blinky
       [ ] timer tick function to measure 10ms-ish type times
       [ ] switch debouncing?
       [ ] light led with bit errors
       [ ] 16 to 8 kHz rate conversion
       [ ] change freedv_api interface to float[]
    */

    /* clear filter memories */

    for(i=0; i<FDMDV_OS_TAPS_16K; i++)
	adc16k[i] = 0.0;
    for(i=0; i<FDMDV_OS_TAPS_8K; i++)
	dac8k[i] = 0.0;
    
    while(1) {

        if (0) {

            /* Transmit -------------------------------------------------------------------------*/

            /* ADC2 is the SM1000 microphone, DAC1 is the modulator signal we send to radio tx */

            if (adc1_read(&adc16k[FDMDV_OS_TAPS_16K], FREEDV_NSAMPLES_16K) == 0) {
                GPIOE->ODR = (1 << 3);

                fdmdv_16_to_8_short(adc8k, &adc16k[FDMDV_OS_TAPS_16K], FREEDV_NSAMPLES);

                freedv_tx(f, &dac8k[FDMDV_OS_TAPS_8K], adc8k);
                //for(i=0; i<FREEDV_NSAMPLES; i++)
                //    dac8k[FDMDV_OS_TAPS_8K+i] = adc8k[i];

                fdmdv_8_to_16_short(dac16k, &dac8k[FDMDV_OS_TAPS_8K], FREEDV_NSAMPLES);              

                dac2_write(dac16k, FREEDV_NSAMPLES_16K);

                //led_ptt(1); led_rt(0); led_err(0);
                GPIOE->ODR &= ~(1 << 3);
            }

        }
        else {
            
            /* Receive --------------------------------------------------------------------------*/

            /* ADC1 is the demod in signal from the radio rx, DAC2 is the SM1000 speaker */

            nin = freedv_nin(f);   
            nout = nin;
            f->total_bit_errors = 0;

            if (adc1_read(&adc16k[FDMDV_OS_TAPS_16K], 2*nin) == 0) {
                GPIOE->ODR = (1 << 3);
                fdmdv_16_to_8_short(adc8k, &adc16k[FDMDV_OS_TAPS_16K], nin);
                nout = freedv_rx(f, &dac8k[FDMDV_OS_TAPS_8K], adc8k);
                //for(i=0; i<FREEDV_NSAMPLES; i++)
                //   dac8k[FDMDV_OS_TAPS_8K+i] = adc8k[i];
                fdmdv_8_to_16_short(dac16k, &dac8k[FDMDV_OS_TAPS_8K], nout);              
                dac2_write(dac16k, 2*nout);
                //led_ptt(0); led_rt(f->fdmdv_stats.sync); led_err(f->total_bit_errors);
                GPIOE->ODR &= ~(1 << 3);
            }

        }
       
    } /* while(1) ... */
}

