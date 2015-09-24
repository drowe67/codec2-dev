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
#include "menu.h"

#define FREEDV_NSAMPLES_16K (2*FREEDV_NSAMPLES)

#define MENU_LED_PERIOD  30000
#define ANNOUNCE_DELAY  300000  /* Supposed to be msec, seems not */
#define HOLD_DELAY      100000
#define MENU_DELAY      100000

#define MAX_MODES  3
#define ANALOG     0
#define DV         1
#define TONE       2

struct switch_t sw_select;  /*!< Switch driver for SELECT button */
struct switch_t sw_back;    /*!< Switch driver for BACK button */
struct switch_t sw_ptt;     /*!< Switch driver for PTT buttons */

unsigned int announceTicker = 0;
unsigned int menuLEDTicker = 0;
unsigned int menuTicker = 0;

/*!
 * User preferences
 */
static struct prefs_t {
    /*! Operating mode */
    uint8_t op_mode;
    /*! Menu volume (attenuation) */
    uint8_t menu_vol;
    /* TODO: more to come */
} prefs;

struct tone_gen_t tone_gen;
struct sfx_player_t sfx_player;
struct morse_player_t morse_player;

void SysTick_Handler(void);

/*! Menu item root */
static const struct menu_item_t menu_root;

#define MENU_EVT_NEXT   0x10    /*!< Increment the current item */
#define MENU_EVT_PREV   0x11    /*!< Decrement the current item */
#define MENU_EVT_SELECT 0x20    /*!< Select current item */
#define MENU_EVT_BACK   0x21    /*!< Go back one level */
#define MENU_EVT_EXIT   0x30    /*!< Exit menu */

int main(void) {
    struct freedv *f;
    int            nin, nout, i;
    int            n_samples, n_samples_16k;

    /* Menu data */
    struct menu_t   menu;

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

    /* Clear out switch states */
    memset(&sw_select, 0, sizeof(sw_select));
    memset(&sw_back, 0, sizeof(sw_back));
    memset(&sw_ptt, 0, sizeof(sw_ptt));

    /* Clear out menu state */
    memset(&menu, 0, sizeof(menu));

    morse_player.freq = 800;
    morse_player.dit_time = 60;    /* 20 WPM */
    morse_player.msg = NULL;

    tone_reset(&tone_gen, 0, 0);

    /* play a start-up tune. */
    sfx_play(&sfx_player, sound_startup);

    prefs.op_mode = ANALOG;
    prefs.menu_vol = 7;

    while(1) {
        /* Read switch states */
        switch_update(&sw_select,   (!switch_select()) ? 1 : 0);
        switch_update(&sw_back,     (!switch_back()) ? 1 : 0);
        switch_update(&sw_ptt,      (switch_ptt() ||
                                        (!ext_ptt())) ? 1 : 0);

        /* Process menu */
        if (!menuTicker) {
            if (menu.stack_depth > 0) {
                /* We are in a menu */
                static uint8_t press_ack = 0;

                if (press_ack == 1) {
                    if ((sw_select.state == SW_STEADY)
                            && (!sw_select.sw))
                        press_ack = 0;
                } else if (press_ack == 2) {
                    if ((sw_back.state == SW_STEADY)
                            && (!sw_back.sw))
                        press_ack = 0;
                } else {
                    if (switch_pressed(&sw_select) > HOLD_DELAY) {
                        menu_exec(&menu, MENU_EVT_SELECT);
                        press_ack = 1;
                        menuTicker = MENU_DELAY;
                    } else if (switch_pressed(&sw_back) > HOLD_DELAY) {
                        menu_exec(&menu, MENU_EVT_BACK);
                        press_ack = 2;
                        menuTicker = MENU_DELAY;
                    } else if (switch_released(&sw_select)) {
                        menu_exec(&menu, MENU_EVT_NEXT);
                        menuTicker = MENU_DELAY;
                    } else if (switch_released(&sw_back)) {
                        menu_exec(&menu, MENU_EVT_PREV);
                        menuTicker = MENU_DELAY;
                    } else if (switch_released(&sw_ptt)) {
                        while(menu.stack_depth > 0)
                            menu_exec(&menu, MENU_EVT_EXIT);
                        sfx_play(&sfx_player, sound_returned);
                    }

                    /* If exited, put the LED back */
                    if (!menu.stack_depth) {
                        menuLEDTicker = 0;
                        menuTicker = 0;
                        led_pwr(LED_ON);
                        morse_play(&morse_player, NULL);
                    }
                }
            } else {
                uint8_t mode_changed = 0;

                if (switch_pressed(&sw_select) > HOLD_DELAY) {
                    /* Enter the menu */
                    menu_enter(&menu, &menu_root);
                    menuTicker = MENU_DELAY;
                } else if (switch_released(&sw_select)) {
                    /* Shortcut: change current mode */
                    prefs.op_mode = (prefs.op_mode + 1) % MAX_MODES;
                    mode_changed = 1;
                } else if (switch_released(&sw_back)) {
                    /* Shortcut: change current mode */
                    prefs.op_mode = (prefs.op_mode - 1) % MAX_MODES;
                    mode_changed = 1;
                }

                if (mode_changed) {
                    /* Announce the new mode */
                    if (prefs.op_mode == ANALOG)
                        morse_play(&morse_player, "ANA");
                    else if (prefs.op_mode == DV)
                        morse_play(&morse_player, "1600");
                    else if (prefs.op_mode == TONE)
                        morse_play(&morse_player, "TONE");
                    sfx_play(&sfx_player, sound_click);
                }
            }
        }

        /* Acknowledge switch events */
        switch_ack(&sw_select);
        switch_ack(&sw_back);
        switch_ack(&sw_ptt);

        if (sfx_player.note) {
            int samples = 0;
            int sz_free = dac2_free();
            if (sz_free > n_samples_16k)
                sz_free = n_samples_16k;
            for (i=0; i < sz_free; i++) {
                dac16k[i] = sfx_next(&sfx_player) >> prefs.menu_vol;
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
                dac16k[i] = morse_next(&morse_player) >> prefs.menu_vol;
                samples++;
                if (!morse_player.msg)
                    break;
            }
            dac2_write(dac16k, samples);
        }
        else if (menu.stack_depth > 0) {
            if (!menuLEDTicker) {
                led_pwr(LED_INV);
                menuLEDTicker = MENU_LED_PERIOD;
            }
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

                if (prefs.op_mode == ANALOG) {
                    for(i=0; i<n_samples; i++)
                        dac8k[FDMDV_OS_TAPS_8K+i] = adc8k[i];
                    fdmdv_8_to_16_short(dac16k, &dac8k[FDMDV_OS_TAPS_8K], n_samples);
                    dac1_write(dac16k, n_samples_16k);
                }
                if (prefs.op_mode == DV) {
                    freedv_tx(f, &dac8k[FDMDV_OS_TAPS_8K], adc8k);
                    for(i=0; i<n_samples; i++)
                        dac8k[FDMDV_OS_TAPS_8K+i] *= 0.398; /* 8dB back off from peak */
                    fdmdv_8_to_16_short(dac16k, &dac8k[FDMDV_OS_TAPS_8K], n_samples);
                    dac1_write(dac16k, n_samples_16k);
                }
                if (prefs.op_mode == TONE) {
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

            if (prefs.op_mode == ANALOG) {

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
    switch_tick(&sw_select);
    switch_tick(&sw_back);
    switch_tick(&sw_ptt);
    if (menuTicker > 0) {
        menuTicker--;
    }
    if (menuLEDTicker > 0) {
        menuLEDTicker--;
    }
    if (announceTicker > 0) {
        announceTicker--;
    }
}

/* ---------------------------- Menu data --------------------------- */

/*!
 * Default handler for menu callback.
 */
static void menu_default_cb(struct menu_t* const menu, uint32_t event)
{
    /* Get the current menu item */
    const struct menu_item_t* item = menu_item(menu, 0);
    uint8_t announce = 0;

    switch(event) {
        case MENU_EVT_ENTERED:
            sfx_play(&sfx_player, sound_startup);
            /* Choose first item */
            menu->current = 0;
        case MENU_EVT_RETURNED:
            announce = 1;
            break;
        case MENU_EVT_NEXT:
            sfx_play(&sfx_player, sound_click);
            menu->current = (menu->current + 1) % item->num_children;
            announce = 1;
            break;
        case MENU_EVT_PREV:
            sfx_play(&sfx_player, sound_click);
            menu->current = (menu->current - 1) % item->num_children;
            announce = 1;
            break;
        case MENU_EVT_SELECT:
            /* Enter the sub-menu */
            menu_enter(menu, item->children[menu->current]);
            break;
        case MENU_EVT_BACK:
            /* Exit the menu */
            sfx_play(&sfx_player, sound_returned);
        case MENU_EVT_EXIT:
            menu_leave(menu);
            break;
        default:
            break;
    }

    if (announce) {
        /* Announce the label of the selected child */
        morse_play(&morse_player,
                item->children[menu->current]->label);
    }
}

/* Root menu item forward declarations */
static const struct menu_item_t const* menu_root_children[];
/* Root item definition */
static const struct menu_item_t menu_root = {
    .label          = "MENU",
    .event_cb       = menu_default_cb,
    .children       = menu_root_children,
    .num_children   = 2,
};

/* Child declarations */
static const struct menu_item_t menu_op_mode;
static const struct menu_item_t menu_ui;
static const struct menu_item_t const* menu_root_children[] = {
    &menu_op_mode,
    &menu_ui,
};


/* Operation Mode menu forward declarations */
static void menu_op_mode_cb(struct menu_t* const menu, uint32_t event);
static struct menu_item_t const* menu_op_mode_children[];
/* Operation mode menu */
static const struct menu_item_t menu_op_mode = {
    .label          = "MODE",
    .event_cb       = menu_op_mode_cb,
    .children       = menu_op_mode_children,
    .num_children   = 3,
};
/* Children */
static const struct menu_item_t menu_op_mode_analog = {
    .label          = "ANA",
    .event_cb       = NULL,
    .children       = NULL,
    .num_children   = 0,
    .data           = {
        .ui         = ANALOG,
    },
};
static const struct menu_item_t menu_op_mode_dv16k = {
    .label          = "1600",
    .event_cb       = NULL,
    .children       = NULL,
    .num_children   = 0,
    .data           = {
        .ui         = DV,
    },
};
/* static const struct menu_item_t menu_op_mode_dv700b
    .label          = "700",
    .event_cb       = NULL,
    .children       = NULL,
    .num_children   = 0,
    .data           = {
        .ui         = DV,
    },
};*/
static const struct menu_item_t menu_op_mode_tone = {
    .label          = "TONE",
    .event_cb       = NULL,
    .children       = NULL,
    .num_children   = 0,
    .data           = {
        .ui         = TONE,
    },
};
static struct menu_item_t const* menu_op_mode_children[] = {
    &menu_op_mode_analog,
    &menu_op_mode_dv16k,
    /* &menu_op_mode_dv700b, */
    &menu_op_mode_tone,
};
/* Callback function */
static void menu_op_mode_cb(struct menu_t* const menu, uint32_t event)
{
    const struct menu_item_t* item = menu_item(menu, 0);
    uint8_t announce = 0;

    switch(event) {
        case MENU_EVT_ENTERED:
            sfx_play(&sfx_player, sound_startup);
            /* Choose current item */
            switch(prefs.op_mode) {
                case DV:
                    menu->current = 1;
                    break;
                case TONE:
                    menu->current = 2;
                    break;
                default:
                    menu->current = 0;
            }
        case MENU_EVT_RETURNED:
            /* Shouldn't happen, but we handle it anyway */
            announce = 1;
            break;
        case MENU_EVT_NEXT:
            sfx_play(&sfx_player, sound_click);
            menu->current = (menu->current + 1) % item->num_children;
            announce = 1;
            break;
        case MENU_EVT_PREV:
            sfx_play(&sfx_player, sound_click);
            menu->current = (menu->current - 1) % item->num_children;
            announce = 1;
            break;
        case MENU_EVT_SELECT:
            /* Choose the selected mode */
            prefs.op_mode = item->children[menu->current]->data.ui;
            /* Play the "selected" tune and return. */
            sfx_play(&sfx_player, sound_startup);
            menu_leave(menu);
            break;
        case MENU_EVT_BACK:
            /* Exit the menu */
            sfx_play(&sfx_player, sound_returned);
        case MENU_EVT_EXIT:
            menu_leave(menu);
            break;
        default:
            break;
    }

    if (announce) {
        /* Announce the label of the selected child */
        morse_play(&morse_player,
                item->children[menu->current]->label);
    }
}


/* UI menu forward declarations */
static struct menu_item_t const* menu_ui_children[];
/* Operation mode menu */
static const struct menu_item_t menu_ui = {
    .label          = "UI",
    .event_cb       = menu_default_cb,
    .children       = menu_ui_children,
    .num_children   = 3,
};
/* Children */
static const struct menu_item_t menu_ui_freq;
static const struct menu_item_t menu_ui_speed;
static const struct menu_item_t menu_ui_vol;
static struct menu_item_t const* menu_ui_children[] = {
    &menu_ui_freq,
    &menu_ui_speed,
    &menu_ui_vol,
};

/* UI Frequency menu forward declarations */
static void menu_ui_freq_cb(struct menu_t* const menu, uint32_t event);
/* UI Frequency menu */
static const struct menu_item_t menu_ui_freq = {
    .label          = "FREQ",
    .event_cb       = menu_ui_freq_cb,
    .children       = NULL,
    .num_children   = 0,
};
/* Callback function */
static void menu_ui_freq_cb(struct menu_t* const menu, uint32_t event)
{
    uint8_t announce = 0;

    switch(event) {
        case MENU_EVT_ENTERED:
            sfx_play(&sfx_player, sound_startup);
            /* Get the current frequency */
            menu->current = morse_player.freq;
        case MENU_EVT_RETURNED:
            /* Shouldn't happen, but we handle it anyway */
            announce = 1;
            break;
        case MENU_EVT_NEXT:
            sfx_play(&sfx_player, sound_click);
            /* Adjust the frequency up by 50 Hz */
            if (morse_player.freq < 2000)
                morse_player.freq += 50;
            announce = 1;
            break;
        case MENU_EVT_PREV:
            sfx_play(&sfx_player, sound_click);
            if (morse_player.freq > 50)
                morse_player.freq -= 50;
            announce = 1;
            break;
        case MENU_EVT_SELECT:
            /* Play the "selected" tune and return. */
            sfx_play(&sfx_player, sound_startup);
            menu_leave(menu);
            break;
        case MENU_EVT_BACK:
            /* Restore the mode and exit the menu */
            sfx_play(&sfx_player, sound_returned);
        case MENU_EVT_EXIT:
            morse_player.freq = menu->current;
            menu_leave(menu);
            break;
        default:
            break;
    }

    if (announce) {
        /* Render the text, thankfully we don't need re-entrancy */
        static char freq[5];
        snprintf(freq, 4, "%d", morse_player.freq);
        /* Announce the frequency */
        morse_play(&morse_player, freq);
    }
};

/* UI Speed menu forward declarations */
static void menu_ui_speed_cb(struct menu_t* const menu, uint32_t event);
/* UI Speed menu */
static const struct menu_item_t menu_ui_speed = {
    .label          = "WPM",
    .event_cb       = menu_ui_speed_cb,
    .children       = NULL,
    .num_children   = 0,
};
/* Callback function */
static void menu_ui_speed_cb(struct menu_t* const menu, uint32_t event)
{
    uint8_t announce = 0;

    /* Get the current WPM */
    uint16_t curr_wpm = 1200 / morse_player.dit_time;

    switch(event) {
        case MENU_EVT_ENTERED:
            sfx_play(&sfx_player, sound_startup);
            /* Get the current frequency */
            menu->current = morse_player.dit_time;
        case MENU_EVT_RETURNED:
            /* Shouldn't happen, but we handle it anyway */
            announce = 1;
            break;
        case MENU_EVT_NEXT:
            sfx_play(&sfx_player, sound_click);
            /* Increment WPM by 5 */
            if (curr_wpm < 60)
                curr_wpm += 5;
            announce = 1;
            break;
        case MENU_EVT_PREV:
            sfx_play(&sfx_player, sound_click);
            if (curr_wpm > 5)
                curr_wpm -= 5;
            announce = 1;
            break;
        case MENU_EVT_SELECT:
            /* Play the "selected" tune and return. */
            sfx_play(&sfx_player, sound_startup);
            menu_leave(menu);
            break;
        case MENU_EVT_BACK:
            /* Restore the mode and exit the menu */
            sfx_play(&sfx_player, sound_returned);
        case MENU_EVT_EXIT:
            morse_player.dit_time = menu->current;
            menu_leave(menu);
            break;
        default:
            break;
    }

    if (announce) {
        /* Render the text, thankfully we don't need re-entrancy */
        static char wpm[5];
        snprintf(wpm, 4, "%d", curr_wpm);
        /* Set the new parameter */
        morse_player.dit_time = 1200 / curr_wpm;
        /* Announce the words per minute */
        morse_play(&morse_player, wpm);
    }
};

/* UI volume menu forward declarations */
static void menu_ui_vol_cb(struct menu_t* const menu, uint32_t event);
/* UI volume menu */
static const struct menu_item_t menu_ui_vol = {
    .label          = "VOL",
    .event_cb       = menu_ui_vol_cb,
    .children       = NULL,
    .num_children   = 0,
};
/* Callback function */
static void menu_ui_vol_cb(struct menu_t* const menu, uint32_t event)
{
    uint8_t announce = 0;

    switch(event) {
        case MENU_EVT_ENTERED:
            sfx_play(&sfx_player, sound_startup);
            /* Get the current volume */
            menu->current = prefs.menu_vol;
        case MENU_EVT_RETURNED:
            /* Shouldn't happen, but we handle it anyway */
            announce = 1;
            break;
        case MENU_EVT_NEXT:
            sfx_play(&sfx_player, sound_click);
            if (prefs.menu_vol > 0)
                prefs.menu_vol--;
            announce = 1;
            break;
        case MENU_EVT_PREV:
            sfx_play(&sfx_player, sound_click);
            if (prefs.menu_vol < 14)
                prefs.menu_vol++;
            announce = 1;
            break;
        case MENU_EVT_SELECT:
            /* Play the "selected" tune and return. */
            sfx_play(&sfx_player, sound_startup);
            menu_leave(menu);
            break;
        case MENU_EVT_BACK:
            /* Restore the mode and exit the menu */
            sfx_play(&sfx_player, sound_returned);
        case MENU_EVT_EXIT:
            morse_player.dit_time = menu->current;
            menu_leave(menu);
            break;
        default:
            break;
    }

    if (announce) {
        /* Render the text, thankfully we don't need re-entrancy */
        static char vol[3];
        snprintf(vol, 2, "%d", 15 - prefs.menu_vol);
        /* Announce the volume level */
        morse_play(&morse_player, vol);
    }
};
