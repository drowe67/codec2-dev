/*---------------------------------------------------------------------------*\

  FILE........: timer_ut.c
  AUTHOR......: David Rowe
  DATE CREATED: 3 Jan 2014

  Unit test STM32F4 timer hardware.

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

#include <assert.h>
#include <stdlib.h>
#include <string.h>

#include "stm32f4xx_gpio.h"
#include "stm32f4xx_rcc.h"

#include "gdb_stdio.h"

#define TIM1_CCR3_ADDRESS    0x4001223C

TIM_TimeBaseInitTypeDef  TIM_TimeBaseStructure;
TIM_OCInitTypeDef  TIM_OCInitStructure;
TIM_BDTRInitTypeDef TIM_BDTRInitStructure;
uint16_t uhTimerPeriod;
uint16_t aSRC_Buffer[3] = {0, 0, 0};

void Timer1Config();
#define FS  3500000

int main(void){
    Timer1Config();
 }

/* DR: TIM_Config configures a couple of I/O pins for PWM output from
   Timer1 Channel 3.  Note I dont think any of this is needed, except
   perhaps to check timer frequency.  Can be removed down the track. */

/**
  * @brief  Configure the TIM1 Pins.
  * @param  None
  * @retval None
  */
static void TIM_Config(void)
{
  GPIO_InitTypeDef GPIO_InitStructure;

  /* GPIOA and GPIOB clock enable */
  RCC_AHB1PeriphClockCmd(RCC_AHB1Periph_GPIOA | RCC_AHB1Periph_GPIOB, ENABLE);

  /* GPIOA Configuration: Channel 3 as alternate function push-pull */
  /* Discovery board pin PA10 */

  GPIO_InitStructure.GPIO_Pin = GPIO_Pin_10 ;
  GPIO_InitStructure.GPIO_Mode = GPIO_Mode_AF;
  GPIO_InitStructure.GPIO_Speed = GPIO_Speed_100MHz;
  GPIO_InitStructure.GPIO_OType = GPIO_OType_PP;
  GPIO_InitStructure.GPIO_PuPd = GPIO_PuPd_UP ;
  GPIO_Init(GPIOA, &GPIO_InitStructure);
  GPIO_PinAFConfig(GPIOA, GPIO_PinSource10, GPIO_AF_TIM1);

  /* GPIOB Configuration: Channel 3N as alternate function push-pull */
  /* Discovery board pin PB15 */

  GPIO_InitStructure.GPIO_Pin = GPIO_Pin_15;
  GPIO_Init(GPIOB, &GPIO_InitStructure);
  GPIO_PinAFConfig(GPIOB, GPIO_PinSource15, GPIO_AF_TIM1);
}

void Timer1Config() {

    /* TIM Configuration */

    TIM_Config();

    /* TIM1 example -------------------------------------------------

       TIM1 input clock (TIM1CLK) is set to 2 * APB2 clock (PCLK2), since APB2
       prescaler is different from 1.
       TIM1CLK = 2 * PCLK2
       PCLK2 = HCLK / 2
       => TIM1CLK = 2 * (HCLK / 2) = HCLK = SystemCoreClock

       TIM1CLK = SystemCoreClock, Prescaler = 0, TIM1 counter clock = SystemCoreClock
       SystemCoreClock is set to 168 MHz for STM32F4xx devices.

       The objective is to configure TIM1 channel 3 to generate complementary PWM
       signal with a frequency equal to F KHz:
       - TIM1_Period = (SystemCoreClock / F) - 1

       The number of this repetitive requests is defined by the TIM1 Repetion counter,
       each 3 Update Requests, the TIM1 Channel 3 Duty Cycle changes to the next new
       value defined by the aSRC_Buffer.

       Note:
       SystemCoreClock variable holds HCLK frequency and is defined in system_stm32f4xx.c file.
       Each time the core clock (HCLK) changes, user had to call SystemCoreClockUpdate()
       function to update SystemCoreClock variable value. Otherwise, any configuration
       based on this variable will be incorrect.
       -----------------------------------------------------------------------------*/

    /* Compute the value to be set in ARR regiter to generate signal frequency at FS Hz */
    uhTimerPeriod = (SystemCoreClock / FS ) - 1;
    /* Compute CCR1 value to generate a duty cycle at 50% */
    aSRC_Buffer[0] = (uint16_t) (((uint32_t) 5 * (uhTimerPeriod - 1)) / 10);
    /* Compute CCR1 value to generate a duty cycle at 37.5% */
    aSRC_Buffer[1] = (uint16_t) (((uint32_t) 375 * (uhTimerPeriod - 1)) / 1000);
    /* Compute CCR1 value to generate a duty cycle at 25% */
    aSRC_Buffer[2] = (uint16_t) (((uint32_t) 25 * (uhTimerPeriod - 1)) / 100);

    /* TIM1 Peripheral Configuration -------------------------------------------*/
    /* TIM1 clock enable */
    RCC_APB2PeriphClockCmd(RCC_APB2Periph_TIM1, ENABLE);

    /* Time Base configuration */

    TIM_DeInit(TIM1);
    TIM_TimeBaseStructure.TIM_Prescaler = 0;
    TIM_TimeBaseStructure.TIM_CounterMode = TIM_CounterMode_Up;
    TIM_TimeBaseStructure.TIM_Period = uhTimerPeriod;
    TIM_TimeBaseStructure.TIM_ClockDivision = 0;
    TIM_TimeBaseStructure.TIM_RepetitionCounter = 0;

    TIM_TimeBaseInit(TIM1, &TIM_TimeBaseStructure);

    /* Channel 3 Configuration in PWM mode */

    /* I think we just ned to enable channel 3 somehow, but without
       (or optionally with) actual ouput to a GPIO pin.  */

    TIM_OCInitStructure.TIM_OCMode = TIM_OCMode_PWM2;
    TIM_OCInitStructure.TIM_OutputState = TIM_OutputState_Enable;
    TIM_OCInitStructure.TIM_OutputNState = TIM_OutputNState_Enable;
    TIM_OCInitStructure.TIM_Pulse = aSRC_Buffer[0];
    TIM_OCInitStructure.TIM_OCPolarity = TIM_OCPolarity_Low;
    TIM_OCInitStructure.TIM_OCNPolarity = TIM_OCNPolarity_Low;
    TIM_OCInitStructure.TIM_OCIdleState = TIM_OCIdleState_Set;
    TIM_OCInitStructure.TIM_OCNIdleState = TIM_OCIdleState_Reset;

    TIM_OC3Init(TIM1, &TIM_OCInitStructure);

    /* Enable preload feature */
    TIM_OC3PreloadConfig(TIM1, TIM_OCPreload_Enable);

    /* Automatic Output enable, Break, dead time and lock configuration*/
    TIM_BDTRInitStructure.TIM_OSSRState = TIM_OSSRState_Enable;
    TIM_BDTRInitStructure.TIM_OSSIState = TIM_OSSIState_Enable;
    //TIM_BDTRInitStructure.TIM_LOCKLevel = TIM_LOCKLevel_1;
    TIM_BDTRInitStructure.TIM_DeadTime = 11;
    //TIM_BDTRInitStructure.TIM_Break = TIM_Break_Enable;
    //TIM_BDTRInitStructure.TIM_BreakPolarity = TIM_BreakPolarity_High;
    TIM_BDTRInitStructure.TIM_AutomaticOutput = TIM_AutomaticOutput_Enable;

    TIM_BDTRConfig(TIM1, &TIM_BDTRInitStructure);

    /* TIM1 counter enable */
    TIM_Cmd(TIM1, ENABLE);

    /* Main Output Enable */
    TIM_CtrlPWMOutputs(TIM1, ENABLE);
}

