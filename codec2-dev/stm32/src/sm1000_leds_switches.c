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

#define LED_PWR       12
#define LED_PTT       13
#define LED_RT        14
#define LED_ERR       15
#define SWITCH_PTT     7
#define SWITCH_SELECT  0
#define SWITCH_BACK    1

#include <stm32f4xx.h>
#include <stm32f4xx_gpio.h>
#include "sm1000_leds_switches.h"

void sm1000_leds_switches_init(void) {
    GPIO_InitTypeDef GPIO_InitStruct;

    RCC_AHB1PeriphClockCmd(RCC_AHB1Periph_GPIOD, ENABLE);

    /* output pins */

    GPIO_InitStruct.GPIO_Pin = LED_PWR | LED_PTT | LED_RT | LED_ERR;
    GPIO_InitStruct.GPIO_Mode = GPIO_Mode_OUT; 		
    GPIO_InitStruct.GPIO_Speed = GPIO_Speed_50MHz; 	
    GPIO_InitStruct.GPIO_OType = GPIO_OType_PP; 	 
    GPIO_InitStruct.GPIO_PuPd = GPIO_PuPd_NOPULL; 	
    GPIO_Init(GPIOD, &GPIO_InitStruct); 		

    /* input pins */

    GPIO_InitStruct.GPIO_Pin = SWITCH_PTT | SWITCH_SELECT | SWITCH_BACK;
    GPIO_InitStruct.GPIO_Mode = GPIO_Mode_IN; 		
    GPIO_InitStruct.GPIO_Speed = GPIO_Speed_50MHz; 
    GPIO_InitStruct.GPIO_PuPd = GPIO_PuPd_NOPULL; /* we have our own external pull ups */	
    GPIO_Init(GPIOD, &GPIO_InitStruct); 		
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

int switch_ptt(void) {
    return GPIOA->IDR & SWITCH_PTT;
}

int switch_select(void) {
    return GPIOA->IDR & SWITCH_SELECT;
}

int switch_back(void) {
    return GPIOA->IDR & SWITCH_BACK;
}
