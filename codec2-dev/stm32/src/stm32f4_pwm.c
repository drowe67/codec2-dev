/*---------------------------------------------------------------------------*\

  FILE........: stm32f4_pwm.c
  AUTHOR......: David Rowe
  DATE CREATED: 26 June 2013

  PWM  driver module for STM32F4.

  TODO:

\*---------------------------------------------------------------------------*/

/*
  Copyright (C) 2013 David Rowe

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
 
#define TIM1_CCR3_ADDRESS    0x4001003C
//#define TIM1_CCR3_ADDRESS    0x4001223C
#define SINE_SAMPLES         32

TIM_TimeBaseInitTypeDef  TIM_TimeBaseStructure;
TIM_OCInitTypeDef  TIM_OCInitStructure;
uint16_t uhTimerPeriod;
uint16_t aSRC_Buffer[SINE_SAMPLES] = {0, 0, 0};

/* 32 sample sine wave which at Fs=16kHz will be 500Hz.  Not sampels
   are 16 bit 2's complement, the DAC driver convertsto 12 bit
   unsigned. */

short aSine[SINE_SAMPLES] = {
    -16,    6384,   12528,  18192,   23200,   27232,   30256,   32128,   32752,   32128,
    30256,   27232,   23152,   18192,   12528,    6384,     -16,   -6416,  -12560,  -18224,
    -23184,  -27264,  -30288,  -32160,  -32768,  -32160,  -30288,  -27264,  -23184,  -18224,
    -12560,   -6416
};

void Timer1Config();

#define FS  16000

int main(void){
    Timer1Config();
    while(1);
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
  DMA_InitTypeDef DMA_InitStructure;
  
  /* GPIOA and GPIOB clock enable */
  RCC_AHB1PeriphClockCmd(RCC_AHB1Periph_GPIOA | RCC_AHB1Periph_GPIOB, ENABLE);

  /* GPIOA Configuration: Channel 3 as alternate function push-pull */

  GPIO_InitStructure.GPIO_Pin = GPIO_Pin_10 ;
  GPIO_InitStructure.GPIO_Mode = GPIO_Mode_AF;
  GPIO_InitStructure.GPIO_Speed = GPIO_Speed_100MHz;
  GPIO_InitStructure.GPIO_OType = GPIO_OType_PP;
  GPIO_InitStructure.GPIO_PuPd = GPIO_PuPd_UP ;
  GPIO_Init(GPIOA, &GPIO_InitStructure); 
  GPIO_PinAFConfig(GPIOA, GPIO_PinSource10, GPIO_AF_TIM1);

  /* GPIOB Configuration: Channel 3N as alternate function push-pull */

  GPIO_InitStructure.GPIO_Pin = GPIO_Pin_15;
  GPIO_Init(GPIOB, &GPIO_InitStructure);
  GPIO_PinAFConfig(GPIOB, GPIO_PinSource15, GPIO_AF_TIM1);

  /* DMA clock enable */
  RCC_AHB1PeriphClockCmd(RCC_AHB1Periph_DMA2 , ENABLE);

  DMA_DeInit(DMA2_Stream6);
  DMA_InitStructure.DMA_Channel = DMA_Channel_6;  
  DMA_InitStructure.DMA_PeripheralBaseAddr = (uint32_t)(TIM1_CCR3_ADDRESS) ;
  DMA_InitStructure.DMA_Memory0BaseAddr = (uint32_t)aSRC_Buffer;
  DMA_InitStructure.DMA_DIR = DMA_DIR_MemoryToPeripheral;
  DMA_InitStructure.DMA_BufferSize = SINE_SAMPLES;
  DMA_InitStructure.DMA_PeripheralInc = DMA_PeripheralInc_Disable;
  DMA_InitStructure.DMA_MemoryInc = DMA_MemoryInc_Enable;
  DMA_InitStructure.DMA_PeripheralDataSize = DMA_PeripheralDataSize_HalfWord;
  DMA_InitStructure.DMA_MemoryDataSize = DMA_PeripheralDataSize_HalfWord;
  DMA_InitStructure.DMA_Mode = DMA_Mode_Circular;
  DMA_InitStructure.DMA_Priority = DMA_Priority_High;
  DMA_InitStructure.DMA_FIFOMode = DMA_FIFOMode_Disable;
  DMA_InitStructure.DMA_FIFOThreshold = DMA_FIFOThreshold_Full;
  DMA_InitStructure.DMA_MemoryBurst = DMA_MemoryBurst_Single;
  DMA_InitStructure.DMA_PeripheralBurst = DMA_PeripheralBurst_Single;

  DMA_Init(DMA2_Stream6, &DMA_InitStructure);
}

void Timer1Config() {
    int i;

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
  
    /* Compute the value to be set in ARR regiter to generate signal frequency at FS */

    uhTimerPeriod = (SystemCoreClock / FS ) - 1;

    /* Compute CCR1 values to generate a duty cycle at 50% */

    for(i=0; i<SINE_SAMPLES; i++) {
        aSRC_Buffer[i] = uhTimerPeriod *((int)aSine[i] + 32768)/(32768*2);
    }

#ifdef OLD
  /* Compute CCR1 value to generate a duty cycle at 50% */
  aSRC_Buffer[0] = (uint16_t) (((uint32_t) 5 * (uhTimerPeriod - 1)) / 10);
  /* Compute CCR1 value to generate a duty cycle at 37.5% */
  aSRC_Buffer[1] = (uint16_t) (((uint32_t) 375 * (uhTimerPeriod - 1)) / 1000);
  /* Compute CCR1 value to generate a duty cycle at 25% */
  aSRC_Buffer[2] = (uint16_t) (((uint32_t) 25 * (uhTimerPeriod - 1)) / 100);
#endif

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
  
    /* TIM1 counter enable */
    TIM_Cmd(TIM1, ENABLE);
  
    /* DMA enable*/
    DMA_Cmd(DMA2_Stream6, ENABLE);
  
    /* TIM1 Update DMA Request enable */
    TIM_DMACmd(TIM1, TIM_DMA_CC3, ENABLE);

    /* Main Output Enable */
    TIM_CtrlPWMOutputs(TIM1, ENABLE);
}
