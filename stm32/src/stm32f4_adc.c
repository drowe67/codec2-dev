/*---------------------------------------------------------------------------*\

  FILE........: stm32f4_adc.c
  AUTHOR......: David Rowe
  DATE CREATED: 4 June 2013

  ADC driver module for STM32F4.

  TODO:
  [X] just get ADC to run at all, prove its sampling something....
  [X] as above with DMA
  [X] half and finished interrupts, ISR
  [X] timer config to drive ADC conversion, measure sample rate and confirm 16kHz
      + larger ADC DMA buffer
      + fifos
      + work out a way to unit test
  [ ] ADC working at same time as DAC
  [X] remove (or make optional) the TIM_Config() code that sends PWM output to pins
  [ ] check comments still valid
  [X] convert to driver
  [ ] way to determine which timers are used so they don't get re-sued
  [ ] way to select different pins/ADCs for multiple channels, multiple channel support
  [ ] access functions for halff/full/overflow to trap any issues
  [ ] should FIFOs be in this drivr or in UTs connected to stdio?  SmartMic will just need
      40ms of buffering

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

#include "stm32f4xx_adc.h"
#include "stm32f4xx_gpio.h"
#include "stm32f4xx_rcc.h"
 
#include "codec2_fifo.h"
#include "gdb_stdio.h"
#include "stm32f4_adc.h"

#define ADC_BUF_SZ   320
#define FIFO_SZ      1000

struct FIFO *DMA2_Stream0_fifo;
unsigned short adc_buf[ADC_BUF_SZ];
int adc_overflow;
int half,full;

#define ADCx_DR_ADDRESS          ((uint32_t)0x4001204C)
#define DMA_CHANNELx             DMA_Channel_0
#define DMA_STREAMx              DMA2_Stream0
#define ADCx                     ADC1

#define TIM1_CCR3_ADDRESS    0x4001223C

TIM_TimeBaseInitTypeDef  TIM_TimeBaseStructure;
TIM_OCInitTypeDef  TIM_OCInitStructure;
uint16_t uhTimerPeriod;
uint16_t aSRC_Buffer[3] = {0, 0, 0};

void Timer1Config();
void adc_configure();

void adc_open(void) {
    DMA2_Stream0_fifo = fifo_create(FIFO_SZ);
    assert(DMA2_Stream0_fifo != NULL);

    Timer1Config();
    adc_configure();
    ADC_SoftwareStartConv(ADC1);
}

/* n signed 16 bit samples in buf[] if return != -1 */

int adc_read(short buf[], int n) {   
    return fifo_read(DMA2_Stream0_fifo, buf, n);
}

void Timer1Config() {

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
  
    /* Compute the value to be set in ARR regiter to generate signal frequency at 16.00 Khz */
    uhTimerPeriod = (SystemCoreClock / 16000 ) - 1;
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
  
    /* TIM1 counter enable */
    TIM_Cmd(TIM1, ENABLE);
  
    /* Main Output Enable */
    TIM_CtrlPWMOutputs(TIM1, ENABLE);
}

void adc_configure(){
    ADC_InitTypeDef  ADC_init_structure; 
    GPIO_InitTypeDef GPIO_initStructre; 
    DMA_InitTypeDef  DMA_InitStructure;
    NVIC_InitTypeDef NVIC_InitStructure;

    // Clock configuration

    RCC_APB2PeriphClockCmd(RCC_APB2Periph_ADC1,ENABLE);
    RCC_AHB1PeriphClockCmd(RCC_AHB1ENR_GPIOCEN,ENABLE);
    RCC_AHB1PeriphClockCmd(RCC_AHB1Periph_DMA2, ENABLE);

    // Analog pin configuration

    GPIO_initStructre.GPIO_Pin = GPIO_Pin_0;        // ADC Channel 10 is connected to PC0
    GPIO_initStructre.GPIO_Mode = GPIO_Mode_AN;     
    GPIO_initStructre.GPIO_PuPd = GPIO_PuPd_NOPULL; 
    GPIO_Init(GPIOC,&GPIO_initStructre);            

    // ADC structure configuration

    ADC_DeInit();
    ADC_init_structure.ADC_DataAlign = ADC_DataAlign_Left;
    ADC_init_structure.ADC_Resolution = ADC_Resolution_12b;
    ADC_init_structure.ADC_ContinuousConvMode = DISABLE; 
    ADC_init_structure.ADC_ExternalTrigConv = ADC_ExternalTrigConv_T1_CC3;
    ADC_init_structure.ADC_ExternalTrigConvEdge = ADC_ExternalTrigConvEdge_Rising;
    ADC_init_structure.ADC_NbrOfConversion = 1;
    ADC_init_structure.ADC_ScanConvMode = DISABLE;
    ADC_Init(ADCx,&ADC_init_structure);

    // Select the channel to be read from

    ADC_RegularChannelConfig(ADCx,ADC_Channel_10,1,ADC_SampleTime_144Cycles);

    /* DMA  configuration **************************************/

    DMA_DeInit(DMA_STREAMx);
    DMA_InitStructure.DMA_Channel = DMA_CHANNELx;  
    DMA_InitStructure.DMA_PeripheralBaseAddr = (uint32_t)ADCx_DR_ADDRESS;
    DMA_InitStructure.DMA_Memory0BaseAddr = (uint32_t)adc_buf;
    DMA_InitStructure.DMA_DIR = DMA_DIR_PeripheralToMemory;
    DMA_InitStructure.DMA_BufferSize = ADC_BUF_SZ;
    DMA_InitStructure.DMA_PeripheralInc = DMA_PeripheralInc_Disable;
    DMA_InitStructure.DMA_MemoryInc = DMA_MemoryInc_Enable;
    DMA_InitStructure.DMA_PeripheralDataSize = DMA_PeripheralDataSize_HalfWord;
    DMA_InitStructure.DMA_MemoryDataSize = DMA_MemoryDataSize_HalfWord;
    DMA_InitStructure.DMA_Mode = DMA_Mode_Circular;
    DMA_InitStructure.DMA_Priority = DMA_Priority_High;
    DMA_InitStructure.DMA_FIFOMode = DMA_FIFOMode_Disable;         
    DMA_InitStructure.DMA_FIFOThreshold = DMA_FIFOThreshold_HalfFull;
    DMA_InitStructure.DMA_MemoryBurst = DMA_MemoryBurst_Single;
    DMA_InitStructure.DMA_PeripheralBurst = DMA_PeripheralBurst_Single;
    DMA_Init(DMA_STREAMx, &DMA_InitStructure);

    /* Enable DMA request after last transfer (Single-ADC mode) */

    ADC_DMARequestAfterLastTransferCmd(ADCx, ENABLE);

    /* Enable ADC1 DMA */

    ADC_DMACmd(ADCx, ENABLE);

    /* DMA2_Stream0 enable */

    DMA_Cmd(DMA_STREAMx, ENABLE);

    /* Enable DMA Half & Complete interrupts */

    DMA_ITConfig(DMA2_Stream0, DMA_IT_TC | DMA_IT_HT, ENABLE);

    /* Enable the DMA Stream IRQ Channel */

    NVIC_InitStructure.NVIC_IRQChannel = DMA2_Stream0_IRQn;
    NVIC_InitStructure.NVIC_IRQChannelPreemptionPriority = 0;
    NVIC_InitStructure.NVIC_IRQChannelSubPriority = 0;
    NVIC_InitStructure.NVIC_IRQChannelCmd = ENABLE;
    NVIC_Init(&NVIC_InitStructure);     

    // Enable ADC conversion

    ADC_Cmd(ADC1,ENABLE);
}

/*
  This function handles DMA Stream interrupt request.
*/

void DMA2_Stream0_IRQHandler(void) {
    int i, sam;
    short signed_buf[ADC_BUF_SZ/2];

    /* Half transfer interrupt */

    if(DMA_GetITStatus(DMA2_Stream0, DMA_IT_HTIF0) != RESET) {
        half++;

        /* convert to signed */

        for(i=0; i<ADC_BUF_SZ/2; i++) {
            sam = (int)adc_buf[i] - 32768;
            //sam = (int)adc_buf[i];
            signed_buf[i] = sam;
        }

       /* write first half to fifo */

        if (fifo_write(DMA2_Stream0_fifo, signed_buf, ADC_BUF_SZ/2) == -1) {
            adc_overflow++;
        }

        /* Clear DMA Stream Transfer Complete interrupt pending bit */

        DMA_ClearITPendingBit(DMA2_Stream0, DMA_IT_HTIF0);  
    }

    /* Transfer complete interrupt */

    if(DMA_GetITStatus(DMA2_Stream0, DMA_IT_TCIF0) != RESET) {
        full++;

        /* convert to signed */

        for(i=0; i<ADC_BUF_SZ/2; i++) {
            sam = (int)adc_buf[ADC_BUF_SZ/2 + i] - 32768;
            //sam = (int)adc_buf[ADC_BUF_SZ/2 + i];
            signed_buf[i] = sam;
        }

        /* write second half to fifo */

        if (fifo_write(DMA2_Stream0_fifo, signed_buf, ADC_BUF_SZ/2) == -1) {
            adc_overflow++;
        }

        /* Clear DMA Stream Transfer Complete interrupt pending bit */

        DMA_ClearITPendingBit(DMA2_Stream0, DMA_IT_TCIF0);  
    }
}

