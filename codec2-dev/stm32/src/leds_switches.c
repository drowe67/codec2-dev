/*---------------------------------------------------------------------------*\

  FILE........: leds_switches.c
  AUTHOR......: David Rowe
  DATE CREATED: 18 July 2014

  Functions for controlling LEDs and reading switches on SM1000.

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

#define LED_PWR 12
#define LED_PTT 13
#define LED_RT  14
#define LED_ERR 15

#include "stm32f4xx_conf.h"
#include "stm32f4xx.h"
#include "leds_switches.h"

void leds_switches_init(void) {
    RCC->AHB1ENR |= RCC_AHB1ENR_GPIODEN; // enable the clock to GPIOD 

    // Set pins as general purpose IOs

    GPIOD->MODER = (2 << LED_PWR) | (2 << LED_PTT) | (2 << LED_RT) | (2 << LED_ERR);           
}

void led_pwr(int state) {
    if (state)
        GPIOD->ODR = (1 << LED_PWR);
    else
        GPIOD->ODR &= ~(1 << LED_PWR);
}

void led_ptt(int state) {
    if (state)
        GPIOD->ODR = (1 << LED_PTT);
    else
        GPIOD->ODR &= ~(1 << LED_PTT);
}

void led_rt(int state) {
    if (state)
        GPIOD->ODR = (1 << LED_RT);
    else
        GPIOD->ODR &= ~(1 << LED_RT);
}

void led_err(int state) {
    if (state)
        GPIOD->ODR = (1 << LED_ERR);
    else
        GPIOD->ODR &= ~(1 << LED_ERR);
}

