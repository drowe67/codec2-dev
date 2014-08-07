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
#include "sm1000_leds_switches.h"

int main(void) {
    struct freedv *f;
    short  buf[FREEDV_NSAMPLES];
    int    nin, nout;

    /* init all the drivers for various peripherals */

    sm1000_leds_switches_init();
    dac_open(4*DAC_BUF_SZ);
    adc_open();
    f = freedv_open(FREEDV_MODE_1600);

    /* LEDs into a known state */

    led_pwr(1); led_ptt(0); led_rt(0); led_err(0);

    /* 
       TODO:
       [ ] UT analog interfaces from file IO
       [ ] UTs for simultaneous tx & rx on analog interfaces
       [ ] measure CPU load of various parts with a blinky
       [ ] detect program assert type errors with a blinky
       [ ] timer tick function to measure 10ms-ish type times
       [ ] switch debouncing?
       [ ] light led with bit errors
    */

    while(1) {

        if(switch_ptt()) {

            /* Transmit -------------------------------------------------------------------------*/

            /* ADC2 is the SM1000 microphone, DAC1 is the modulator signal we send to radio tx */

            if (adc2_read(buf, FREEDV_NSAMPLES) == FREEDV_NSAMPLES) {
                freedv_tx(f, buf, buf);
                dac1_write(buf, FREEDV_NSAMPLES);
                led_ptt(1); led_rt(0); led_err(0);
            }
        }
        else {
            
            /* Receive --------------------------------------------------------------------------*/

            /* ADC1 is the demod in signal from the radio rx, DAC2 is the SM1000 speaker */

            nin = freedv_nin(f);
            f->total_bit_errors = 0;
            
            if (adc1_read(buf, nin) == nin) {
                nout = freedv_rx(f, buf, buf);
                dac2_write(buf, nout);
                led_ptt(0); led_rt(f->fdmdv_stats.sync); led_err(f->total_bit_errors);
            }

        }
        
    } /* while(1) ... */
}

