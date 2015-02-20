/*---------------------------------------------------------------------------*\

  FILE........: sm1000_main.c
  AUTHOR......: David Rowe
  DATE CREATED: August 5 2014

  Main program for SM1000.

  TODO

  [ ] make led blink 1-2-3 times for "mode"

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
#include <stdlib.h>

#define FREEDV_NSAMPLES_16K (2*FREEDV_NSAMPLES)

#define FIFTY_MS  50
#define MAX_MODES  3
#define ANALOG     0
#define DV         1
#define TONE       2

#define SS_IDLE           0
#define SS_DEBOUNCE_DOWN  1
#define SS_WAIT_BUTTON_UP 2
#define SS_DEBOUNCE_UP    3

typedef struct {
    int state;
    int mode;
} SWITCH_STATE;

unsigned int downTicker;

void SysTick_Handler(void);
void iterate_select_state_machine(SWITCH_STATE *ss);

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
    struct freedv *f;
    short          adc16k[FDMDV_OS_TAPS_16K+FREEDV_NSAMPLES_16K];
    short          dac16k[FREEDV_NSAMPLES_16K];
    short          adc8k[FREEDV_NSAMPLES];
    short          dac8k[FDMDV_OS_TAPS_8K+FREEDV_NSAMPLES];
    SWITCH_STATE   ss;
    int            nin, nout, i;

    /* init all the drivers for various peripherals */

    SysTick_Config(SystemCoreClock/168000); /* 1 kHz SysTick */
    sm1000_leds_switches_init();
    dac_open(4*DAC_BUF_SZ);
    adc_open(4*ADC_BUF_SZ);
    f = freedv_open(FREEDV_MODE_1600);

    /* put outputs into a known state */

    led_pwr(1); led_ptt(0); led_rt(0); led_err(0); not_cptt(1);

    /* clear filter memories */

    for(i=0; i<FDMDV_OS_TAPS_16K; i++)
	adc16k[i] = 0.0;
    for(i=0; i<FDMDV_OS_TAPS_8K; i++)
	dac8k[i] = 0.0;

    ss.state = SS_IDLE;
    ss.mode  = ANALOG;

    while(1) {
        
        iterate_select_state_machine(&ss);

        if (switch_ptt() || (ext_ptt() == 0)) {
            
            /* Transmit -------------------------------------------------------------------------*/

            /* ADC2 is the SM1000 microphone, DAC1 is the modulator signal we send to radio tx */

            if (adc2_read(&adc16k[FDMDV_OS_TAPS_16K], FREEDV_NSAMPLES_16K) == 0) {
                GPIOE->ODR = (1 << 3);

                /* clipping indicator */

                led_err(0);
                for (i=0; i<FREEDV_NSAMPLES_16K; i++) {
                    if (abs(adc16k[FDMDV_OS_TAPS_16K+i]) > 28000)
                        led_err(1);
                }

                fdmdv_16_to_8_short(adc8k, &adc16k[FDMDV_OS_TAPS_16K], FREEDV_NSAMPLES);

                if (ss.mode == ANALOG) {
                    for(i=0; i<FREEDV_NSAMPLES; i++)
                        dac8k[FDMDV_OS_TAPS_8K+i] = adc8k[i];
                    fdmdv_8_to_16_short(dac16k, &dac8k[FDMDV_OS_TAPS_8K], FREEDV_NSAMPLES);              
                    dac1_write(dac16k, FREEDV_NSAMPLES_16K);
                }
                if (ss.mode == DV) {
                    freedv_tx(f, &dac8k[FDMDV_OS_TAPS_8K], adc8k);
                    for(i=0; i<FREEDV_NSAMPLES; i++)
                        dac8k[FDMDV_OS_TAPS_8K+i] *= 0.398; /* 8dB back off from peak */
                    fdmdv_8_to_16_short(dac16k, &dac8k[FDMDV_OS_TAPS_8K], FREEDV_NSAMPLES);              
                    dac1_write(dac16k, FREEDV_NSAMPLES_16K);
                }
                if (ss.mode == TONE) {
                    short buf[SINE_SAMPLES];
                    for(i=0; i<FREEDV_NSAMPLES; i++)
                        buf[i] = aSine[i]*0.398; /* 8dB back off from peak */                   
                    while(dac1_write(buf, SINE_SAMPLES) == 0);
                }

                led_ptt(1); led_rt(0); led_err(0); not_cptt(0);
                GPIOE->ODR &= ~(1 << 3);
            }

        }
        else {
            
            /* Receive --------------------------------------------------------------------------*/

            not_cptt(1); led_ptt(0); 

            /* ADC1 is the demod in signal from the radio rx, DAC2 is the SM1000 speaker */

            if (ss.mode == ANALOG) {

                /* force analog bypass when select down */

                if (adc1_read(&adc16k[FDMDV_OS_TAPS_16K], FREEDV_NSAMPLES_16K) == 0) {
                    fdmdv_16_to_8_short(adc8k, &adc16k[FDMDV_OS_TAPS_16K], FREEDV_NSAMPLES);
                    for(i=0; i<FREEDV_NSAMPLES; i++)
                        dac8k[FDMDV_OS_TAPS_8K+i] = adc8k[i];
                    fdmdv_8_to_16_short(dac16k, &dac8k[FDMDV_OS_TAPS_8K], FREEDV_NSAMPLES);              
                    dac2_write(dac16k, FREEDV_NSAMPLES_16K);
                    led_rt(0); led_err(0);
               }
            }
            else {

                /* regular DV mode */

                nin = freedv_nin(f);   
                nout = nin;
                f->total_bit_errors = 0;

                if (adc1_read(&adc16k[FDMDV_OS_TAPS_16K], 2*nin) == 0) {
                    GPIOE->ODR = (1 << 3);
                    fdmdv_16_to_8_short(adc8k, &adc16k[FDMDV_OS_TAPS_16K], nin);
                    nout = freedv_rx(f, &dac8k[FDMDV_OS_TAPS_8K], adc8k);
                    fdmdv_8_to_16_short(dac16k, &dac8k[FDMDV_OS_TAPS_8K], nout);              
                    dac2_write(dac16k, 2*nout);
                    led_rt(f->fdmdv_stats.sync); led_err(f->total_bit_errors);
                    GPIOE->ODR &= ~(1 << 3);
                }
            }

        }
    } /* while(1) ... */
}

/*
 * SysTick Interrupt Handler
 */

void SysTick_Handler(void)
{
    if (downTicker > 0) {
        downTicker--;
    }
}

/* Select button state machine.  Debounces switches and enables cycling
   through ANALOG-DV-TONE modes */

void iterate_select_state_machine(SWITCH_STATE *ss) {
    int next_state;

    next_state = ss->state;
    switch(ss->state) {
        case SS_IDLE:
            if (switch_select() == 0) {
                downTicker = FIFTY_MS;
                next_state = SS_DEBOUNCE_DOWN;
            }
            break;
        case SS_DEBOUNCE_DOWN:
            if (downTicker == 0) {
                ss->mode++;
                if (ss->mode >= MAX_MODES)
                    ss->mode = 0;
                next_state = SS_WAIT_BUTTON_UP;
            }
            break;
        case SS_WAIT_BUTTON_UP:
            if (switch_select() == 1) {
                downTicker = FIFTY_MS;
                next_state = SS_DEBOUNCE_UP;
            }
            break;
        case SS_DEBOUNCE_UP:
            if (downTicker == 0) {
                next_state = SS_IDLE;
            }
            break;
   }
    ss->state = next_state;
}
            
