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

#include "sfx.h"
#include "sounds.h"
#include "morse.h"

#define FREEDV_NSAMPLES_16K (2*FREEDV_NSAMPLES)

#define FIFTY_MS  50
#define ANNOUNCE_DELAY  300000  /* Supposed to be msec, seems not */
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

unsigned int downTicker = 0;
unsigned int announceTicker = 0;

struct tone_gen_t tone_gen;
struct sfx_player_t sfx_player;
struct morse_player_t morse_player;

void SysTick_Handler(void);
void iterate_select_state_machine(SWITCH_STATE *ss);

int main(void) {
    struct freedv *f;
    SWITCH_STATE   ss;
    int            nin, nout, i;
    int            n_samples, n_samples_16k;

    /* init all the drivers for various peripherals */

    SysTick_Config(SystemCoreClock/168000); /* 1 kHz SysTick */
    sm1000_leds_switches_init();
    dac_open(4*DAC_BUF_SZ);
    adc_open(4*ADC_BUF_SZ);
    f = freedv_open(FREEDV_MODE_1600);
    n_samples = freedv_get_n_speech_samples(f);
    n_samples_16k = 2*n_samples;

    short          adc16k[FDMDV_OS_TAPS_16K+n_samples_16k];
    short          dac16k[n_samples_16k];
    short          adc8k[n_samples];
    short          dac8k[FDMDV_OS_TAPS_8K+n_samples];

    /* put outputs into a known state */

    led_pwr(1); led_ptt(0); led_rt(0); led_err(0); not_cptt(1);

    /* clear filter memories */

    for(i=0; i<FDMDV_OS_TAPS_16K; i++)
	adc16k[i] = 0.0;
    for(i=0; i<FDMDV_OS_TAPS_8K; i++)
	dac8k[i] = 0.0;

    morse_player.freq = 800;
    morse_player.dit_time = 60;    /* 20 WPM */
    morse_player.msg = NULL;

    tone_reset(&tone_gen, 0, 0);

    /* play a start-up tune. */
    sfx_play(&sfx_player, sound_startup);

    ss.state = SS_IDLE;
    ss.mode  = ANALOG;

    while(1) {

        iterate_select_state_machine(&ss);
        if (sfx_player.note) {
            int samples = 0;
            int sz_free = dac2_free();
            if (sz_free > n_samples_16k)
                sz_free = n_samples_16k;
            for (i=0; i < sz_free; i++) {
                dac16k[i] = sfx_next(&sfx_player) >> 2; /* -6dB */
                samples++;
                if (!sfx_player.note)
                    break;
            }
            dac2_write(dac16k, samples);
            if (!sfx_player.note && morse_player.msg)
                announceTicker = ANNOUNCE_DELAY;
        }
        else if (!announceTicker && morse_player.msg) {
            int samples = 0;
            int sz_free = dac2_free();
            if (sz_free > n_samples_16k)
                sz_free = n_samples_16k;
            for (i=0; i < sz_free; i++) {
                dac16k[i] = morse_next(&morse_player) >> 2; /* -6dB */
                samples++;
                if (!morse_player.msg)
                    break;
            }
            dac2_write(dac16k, samples);
        }
        else if (switch_ptt() || (ext_ptt() == 0)) {

            /* Transmit -------------------------------------------------------------------------*/

            /* Cancel any announcement if scheduled */
            if (announceTicker && morse_player.msg) {
                announceTicker = 0;
                morse_play(&morse_player, NULL);
            }

            /* ADC2 is the SM1000 microphone, DAC1 is the modulator signal we send to radio tx */

            if (adc2_read(&adc16k[FDMDV_OS_TAPS_16K], n_samples_16k) == 0) {
                GPIOE->ODR = (1 << 3);

                /* clipping indicator */

                led_err(0);
                for (i=0; i<n_samples_16k; i++) {
                    if (abs(adc16k[FDMDV_OS_TAPS_16K+i]) > 28000)
                        led_err(1);
                }

                fdmdv_16_to_8_short(adc8k, &adc16k[FDMDV_OS_TAPS_16K], n_samples);

                if (ss.mode == ANALOG) {
                    for(i=0; i<n_samples; i++)
                        dac8k[FDMDV_OS_TAPS_8K+i] = adc8k[i];
                    fdmdv_8_to_16_short(dac16k, &dac8k[FDMDV_OS_TAPS_8K], n_samples);
                    dac1_write(dac16k, n_samples_16k);
                }
                if (ss.mode == DV) {
                    freedv_tx(f, &dac8k[FDMDV_OS_TAPS_8K], adc8k);
                    for(i=0; i<n_samples; i++)
                        dac8k[FDMDV_OS_TAPS_8K+i] *= 0.398; /* 8dB back off from peak */
                    fdmdv_8_to_16_short(dac16k, &dac8k[FDMDV_OS_TAPS_8K], n_samples);
                    dac1_write(dac16k, n_samples_16k);
                }
                if (ss.mode == TONE) {
                    if (!tone_gen.remain)
                        /*
                         * Somewhat ugly, but UINT16_MAX is effectively
                         * infinite.
                         */
                        tone_reset(&tone_gen, 500, UINT16_MAX);
                    int len = dac1_free();
                    if (len > n_samples_16k)
                        len = n_samples_16k;
                    for(i=0; i<len; i++)
                        /* 8dB back off from peak */
                        dac16k[i] = tone_next(&tone_gen)*0.398;
                    dac1_write(dac16k, len);
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

                if (adc1_read(&adc16k[FDMDV_OS_TAPS_16K], n_samples_16k) == 0) {
                    fdmdv_16_to_8_short(adc8k, &adc16k[FDMDV_OS_TAPS_16K], n_samples);
                    for(i=0; i<n_samples; i++)
                        dac8k[FDMDV_OS_TAPS_8K+i] = adc8k[i];
                    fdmdv_8_to_16_short(dac16k, &dac8k[FDMDV_OS_TAPS_8K], n_samples);
                    dac2_write(dac16k, n_samples_16k);
                    led_rt(0); led_err(0);
               }
            }
            else {

                /* regular DV mode */

                nin = freedv_nin(f);
                nout = nin;
		freedv_set_total_bit_errors(f, 0);
                if (adc1_read(&adc16k[FDMDV_OS_TAPS_16K], 2*nin) == 0) {
                    GPIOE->ODR = (1 << 3);
                    fdmdv_16_to_8_short(adc8k, &adc16k[FDMDV_OS_TAPS_16K], nin);
                    nout = freedv_rx(f, &dac8k[FDMDV_OS_TAPS_8K], adc8k);
                    fdmdv_8_to_16_short(dac16k, &dac8k[FDMDV_OS_TAPS_8K], nout);
                    dac2_write(dac16k, 2*nout);
                    led_rt(freedv_get_sync(f)); led_err(freedv_get_total_bit_errors(f));
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
    if (announceTicker > 0) {
        announceTicker--;
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
                sfx_play(&sfx_player, sound_click);
                if (ss->mode >= MAX_MODES)
                    ss->mode = 0;
                if (ss->mode == ANALOG)
                    morse_play(&morse_player, "ANALOG");
                else if (ss->mode == DV)
                    morse_play(&morse_player, "DV");
                else if (ss->mode == TONE)
                    morse_play(&morse_player, "TONE");
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

