/*---------------------------------------------------------------------------*\

  FILE........: sm1000_leds_switches.c
  AUTHOR......: David Rowe
  DATE CREATED: 18 July 2014

  Functions for controlling LEDs and reading switches on the SM1000.

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

#define RIG_PTT        GPIO_Pin_10
#define LED_PWR        GPIO_Pin_12
#define LED_PTT        GPIO_Pin_13
#define LED_SYNC       GPIO_Pin_14
#define LED_ERR        GPIO_Pin_15
#define SWITCH_PTT     GPIO_Pin_7
#define SWITCH_SELECT  GPIO_Pin_0
#define SWITCH_BACK    GPIO_Pin_1
#define EXT_PTT        GPIO_Pin_8

#include <stm32f4xx.h>
#include <stm32f4xx_gpio.h>

#include "sm1000_leds_switches.h"

void sm1000_leds_switches_init(void) {
    GPIO_InitTypeDef GPIO_InitStruct;

    RCC_AHB1PeriphClockCmd(RCC_AHB1Periph_GPIOE, ENABLE);

    /* Debug lights */

    GPIO_InitStruct.GPIO_Pin = GPIO_Pin_0 | GPIO_Pin_1 | GPIO_Pin_2 | GPIO_Pin_3;
    GPIO_InitStruct.GPIO_Mode = GPIO_Mode_OUT;
    GPIO_InitStruct.GPIO_Speed = GPIO_Speed_50MHz;
    GPIO_InitStruct.GPIO_OType = GPIO_OType_PP;
    GPIO_InitStruct.GPIO_PuPd = GPIO_PuPd_NOPULL;

    GPIO_Init(GPIOE, &GPIO_InitStruct);

    RCC_AHB1PeriphClockCmd(RCC_AHB1Periph_GPIOD, ENABLE);

    /* output pins */

    GPIO_InitStruct.GPIO_Pin = LED_PWR | LED_PTT | LED_SYNC | LED_ERR | RIG_PTT;
    GPIO_InitStruct.GPIO_Mode = GPIO_Mode_OUT;
    GPIO_InitStruct.GPIO_Speed = GPIO_Speed_2MHz;
    GPIO_InitStruct.GPIO_OType = GPIO_OType_PP;
    GPIO_InitStruct.GPIO_PuPd = GPIO_PuPd_NOPULL;

    GPIO_Init(GPIOD, &GPIO_InitStruct);

    /* input pins */

    GPIO_InitStruct.GPIO_Pin = SWITCH_PTT | SWITCH_SELECT | SWITCH_BACK;
    GPIO_InitStruct.GPIO_Mode = GPIO_Mode_IN;
    GPIO_InitStruct.GPIO_Speed = GPIO_Speed_2MHz;
    GPIO_InitStruct.GPIO_PuPd = GPIO_PuPd_NOPULL; /* we have our own external pull ups */

    GPIO_Init(GPIOD, &GPIO_InitStruct);

    GPIO_InitStruct.GPIO_Pin = EXT_PTT;
    GPIO_InitStruct.GPIO_Mode = GPIO_Mode_IN;
    GPIO_InitStruct.GPIO_Speed = GPIO_Speed_2MHz;
    GPIO_InitStruct.GPIO_PuPd = GPIO_PuPd_UP;     /* use internal pull up */

    GPIO_Init(GPIOD, &GPIO_InitStruct);
}

/* LED output */

void led_pwr(int state) {
    if (state > 0)
        GPIOD->ODR |= (1 << 12); /*  1 == on */
    else if (state < 0)
        GPIOD->ODR ^= (1 << 12); /* -1 == toggle */
    else
        GPIOD->ODR &= ~(1 << 12);/*  0 == off */
}

void led_ptt(int state) {
    if (state > 0)
        GPIOD->ODR |= (1 << 13); /*  1 == on */
    else if (state < 0)
        GPIOD->ODR ^= (1 << 13); /* -1 == toggle */
    else
        GPIOD->ODR &= ~(1 << 13);/*  0 == off */
}

void led_sync(int state) {
    if (state > 0)
        GPIOD->ODR |= (1 << 14); /*  1 == on */
    else if (state < 0)
        GPIOD->ODR ^= (1 << 14); /* -1 == toggle */
    else
        GPIOD->ODR &= ~(1 << 14);/*  0 == off */
}

void led_err(int state) {
    if (state > 0)
        GPIOD->ODR |= (1 << 15); /*  1 == on */
    else if (state < 0)
        GPIOD->ODR ^= (1 << 15); /* -1 == toggle */
    else
        GPIOD->ODR &= ~(1 << 15);/*  0 == off */
}

void led_debug0(int state) {
    if (state > 0)
        GPIOE->ODR |= (1 << 0);  /*  1 == on */
    else if (state < 0)
        GPIOE->ODR ^= (1 << 0);  /* -1 == toggle */
    else
        GPIOE->ODR &= ~(1 << 0); /*  0 == off */
}

void led_debug1(int state) {
    if (state > 0)
        GPIOE->ODR |= (1 << 1);  /*  1 == on */
    else if (state < 0)
        GPIOE->ODR ^= (1 << 1);  /* -1 == toggle */
    else
        GPIOE->ODR &= ~(1 << 1); /*  0 == off */
}

void led_debug2(int state) {
    if (state > 0)
        GPIOE->ODR |= (1 << 2);  /*  1 == on */
    else if (state < 0)
        GPIOE->ODR ^= (1 << 2);  /* -1 == toggle */
    else
        GPIOE->ODR &= ~(1 << 2); /*  0 == off */
}

void led_debug3(int state) {
    if (state > 0)
        GPIOE->ODR |= (1 << 3);  /*  1 == on */
    else if (state < 0)
        GPIOE->ODR ^= (1 << 3);  /* -1 == toggle */
    else
        GPIOE->ODR &= ~(1 << 3); /*  0 == off */
}

/* Rig Output */

void rig_ptt(int state) {
    /*
     * Flip state for programmer
     */
    if (state)
        GPIOD->ODR &= ~(1 << 10); /* output low for state == high */
    else
        GPIOD->ODR |= (1 << 10);  /* output high for state == low */
}

/* Switch Inputs */

/*
 * Pressing PTT Switch grounds pin
 */
int switch_ptt(void) {
    return GPIOD->IDR & (1 << 7);
}

int switch_select(void) {
    return GPIOD->IDR & (1 << 0);
}

int switch_back(void) {
    return GPIOD->IDR & (1 << 1);
}

/*
 * External PTT Switch grounds pin
 */
int ext_ptt(void) {
    /*
     * Flip state for programmer
     */
    if (GPIOD->IDR & (1 << 8))
        return 0;
    else
        return 1;
}

/*
  FUNCTION: ColorfulRingOfDeath()
  AUTHOR..: xenovacivus

  Colourful ring of death, blink LEDs like crazy forever if something
  really nasty happens.  Adapted from USB Virtual COM Port (VCP)
  module adapted from code I found here:

    https://github.com/xenovacivus/STM32DiscoveryVCP

  Call this to indicate a failure.  Blinks the STM32F4 discovery LEDs
  in sequence.  At 168Mhz, the blinking will be very fast - about 5
  Hz.  Keep that in mind when debugging, knowing the clock speed
  might help with debugging.
*/

int mycode; /* examine this with debugger if it dies */

void ColorfulRingOfDeath(int code) {
    mycode = code;
    uint16_t ring = 1;

    while (1) {
        uint32_t count = 0;

        while (count++ < 5000000)
            ;

        GPIOD->BSRRH = (ring << 12);
        ring = ring << 1;

        if (ring >= 1<<4) {
            ring = 1;
        }

        GPIOD->BSRRL = (ring << 12);
    }
}

void HardFault_Handler(void) { ColorfulRingOfDeath(1); }
void MemManage_Handler(void) { ColorfulRingOfDeath(2); }
void BusFault_Handler(void)  { ColorfulRingOfDeath(3); }
void UsageFault_Handler(void){ ColorfulRingOfDeath(4); }


void switch_tick(struct switch_t* const sw)
{
    if (sw->sw != sw->raw) {
        /* State transition, reset timer */
        if (sw->state == SW_STEADY)
            sw->last = sw->sw;
        sw->state = SW_DEBOUNCE;
        sw->timer = DEBOUNCE_DELAY;
        sw->sw = sw->raw;
    } else if (sw->state == SW_DEBOUNCE) {
        if (sw->timer > 0) {
            /* Steady so far, keep waiting */
            sw->timer--;
        } else {
            /* Steady state reached */
            sw->state = SW_STEADY;
        }
    } else if (sw->sw) {
        /* Hold state.  Yes this will wrap, but who cares? */
        sw->timer++;
    }
}

void switch_update(struct switch_t* const sw, uint8_t state)
{
    sw->raw = state;

    if (sw->raw == sw->sw)
        return;

    if (sw->state == SW_STEADY)
        sw->last = sw->sw;

    sw->timer = DEBOUNCE_DELAY;
    sw->sw = sw->raw;
    sw->state = SW_DEBOUNCE;
}

uint32_t switch_pressed(const struct switch_t* const sw)
{
    if ((sw->state == SW_STEADY) && sw->sw)
        return sw->timer;

    return 0;
}

int switch_released(const struct switch_t* const sw)
{
    if (sw->state != SW_STEADY)
        return 0;

    if (!sw->last)
        return 0;

    if (sw->sw)
        return 0;

    return 1;
}

void switch_ack(struct switch_t* const sw)
{
    if (sw->state == SW_STEADY)
        sw->last = sw->sw;
}
