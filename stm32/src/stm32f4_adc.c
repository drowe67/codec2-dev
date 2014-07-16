/*---------------------------------------------------------------------------*\

  FILE........: stm32f4_adc.c
  AUTHOR......: David Rowe
  DATE CREATED: 4 June 2013

  ADC driver module for STM32F4.

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
#define FIFO_SZ      4*ADC_BUF_SZ

static struct FIFO *adc1_fifo;
static struct FIFO *adc2_fifo;
static unsigned short adc_buf[ADC_BUF_SZ];
static int adc_overflow;
static int half,full;

#define ADCx_DR_ADDRESS          ((uint32_t)0x4001204C)
#define DMA_CHANNELx             DMA_Channel_0
#define DMA_STREAMx              DMA2_Stream0

#define TIM1_CCR3_ADDRESS    0x4001223C

static void Timer1Config();
static void adc_configure();

void adc_open(void) {
    adc1_fifo = fifo_create(FIFO_SZ);
    adc2_fifo = fifo_create(FIFO_SZ);
    assert(adc1_fifo != NULL);
    assert(adc2_fifo != NULL);

    Timer1Config();
    adc_configure();
    ADC_SoftwareStartConv(ADC1);
}

/* n signed 16 bit samples in buf[] if return != -1 */

int adc1_read(short buf[], int n) {   
    return fifo_read(adc1_fifo, buf, n);
}

int adc2_read(short buf[], int n) {   
    return fifo_read(adc2_fifo, buf, n);
}

void Timer1Config() {
    TIM_TimeBaseInitTypeDef  TIM_TimeBaseStructure;
    uint16_t                 uhTimerPeriod;

    /* TIM1 example -------------------------------------------------
  
       TIM1 input clock (TIM1CLK) is set to 2 * APB2 clock (PCLK2), since APB2 
       prescaler is different from 1.   
       TIM1CLK = 2 * PCLK2  
       PCLK2 = HCLK / 2 
       => TIM1CLK = 2 * (HCLK / 2) = HCLK = SystemCoreClock
  
       TIM1CLK = SystemCoreClock, Prescaler = 0, TIM1 counter clock = SystemCoreClock
       SystemCoreClock is set to 168 MHz for STM32F4xx devices.

       The objective is to configure TIM1 channel 3 to generate a
       clock with a frequency equal to F KHz:
       
         TIM1_Period = (SystemCoreClock / F) - 1

       -----------------------------------------------------------------------------*/
  
    /* Compute the value to be set in ARR regiter to generate signal frequency at 16.00 KHz */

    uhTimerPeriod = (SystemCoreClock / 16000 ) - 1;

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

    /* TIM1 counter enable */

    TIM_Cmd(TIM1, ENABLE);
}

void adc_configure() {
    ADC_InitTypeDef  ADC_init_structure; 
    GPIO_InitTypeDef GPIO_initStructre; 
    DMA_InitTypeDef  DMA_InitStructure;
    NVIC_InitTypeDef NVIC_InitStructure;

    // Clock configuration

    RCC_APB2PeriphClockCmd(RCC_APB2Periph_ADC1,  ENABLE);
    RCC_AHB1PeriphClockCmd(RCC_AHB1Periph_GPIOA, ENABLE);
    RCC_AHB1PeriphClockCmd(RCC_AHB1Periph_DMA2,  ENABLE);

    // Analog pin configuration ADC1->PA1, ADC2->PA2 

    GPIO_initStructre.GPIO_Pin = GPIO_Pin_1 | GPIO_Pin_2;        
    GPIO_initStructre.GPIO_Mode = GPIO_Mode_AN;     
    GPIO_initStructre.GPIO_PuPd = GPIO_PuPd_NOPULL; 
    GPIO_Init(GPIOA,&GPIO_initStructre);            

    // ADC structure configuration.  Note we are just using one ADC, which samples
    // two analog channels when triggered by the timer.

    ADC_DeInit();
    ADC_init_structure.ADC_DataAlign = ADC_DataAlign_Left;
    ADC_init_structure.ADC_Resolution = ADC_Resolution_12b;
    ADC_init_structure.ADC_ContinuousConvMode = DISABLE; 
    ADC_init_structure.ADC_ExternalTrigConv = ADC_ExternalTrigConv_T1_CC3;
    ADC_init_structure.ADC_ExternalTrigConvEdge = ADC_ExternalTrigConvEdge_Rising;
    ADC_init_structure.ADC_NbrOfConversion = 2;
    ADC_init_structure.ADC_ScanConvMode = ENABLE;
    ADC_Init(ADC1,&ADC_init_structure);

    // Select the channels to be read from

    ADC_RegularChannelConfig(ADC1,ADC_Channel_1,1,ADC_SampleTime_144Cycles);
    ADC_RegularChannelConfig(ADC1,ADC_Channel_2,2,ADC_SampleTime_144Cycles);

    /* DMA  configuration **************************************/

    DMA_DeInit(DMA2_Stream0);
    DMA_InitStructure.DMA_Channel = DMA_Channel_0;  
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
    DMA_Init(DMA2_Stream0, &DMA_InitStructure);

    /* Enable DMA request after last transfer (Single-ADC mode) */

    ADC_DMARequestAfterLastTransferCmd(ADC1, ENABLE);

    /* Enable ADC1 DMA */

    ADC_DMACmd(ADC1, ENABLE);

    /* DMA2_Stream0 enable */

    DMA_Cmd(DMA2_Stream0, ENABLE);

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
    int i, j, sam1, sam2;
    short signed_buf1[ADC_BUF_SZ/4];
    short signed_buf2[ADC_BUF_SZ/4];

    /* Half transfer interrupt */

    if(DMA_GetITStatus(DMA2_Stream0, DMA_IT_HTIF0) != RESET) {
        half++;

        /* convert to signed */

        for(i=0,j=0; i<ADC_BUF_SZ/4; i++,j+=2) {
            sam1 = (int)adc_buf[j]   - 32768;
            sam2 = (int)adc_buf[j+1] - 32768;
            signed_buf1[i] = sam1;
            signed_buf2[i] = sam2;
        }

        /* write first half to fifo */

        if (fifo_write(adc1_fifo, signed_buf1, ADC_BUF_SZ/4) == -1) {
            adc_overflow++;
        }
        if (fifo_write(adc2_fifo, signed_buf2, ADC_BUF_SZ/4) == -1) {
            adc_overflow++;
        }

        /* Clear DMA Stream Transfer Complete interrupt pending bit */

        DMA_ClearITPendingBit(DMA2_Stream0, DMA_IT_HTIF0);  
    }

    /* Transfer complete interrupt */

    if(DMA_GetITStatus(DMA2_Stream0, DMA_IT_TCIF0) != RESET) {
        full++;

        /* convert to signed */

        for(i=0,j=ADC_BUF_SZ/2; i<ADC_BUF_SZ/4; i++,j+=2) {
            sam1 = (int)adc_buf[j]   - 32768;
            sam2 = (int)adc_buf[j+1] - 32768;
            signed_buf1[i] = sam1;
            signed_buf2[i] = sam2;
        }

        /* write second half to fifo */

        if (fifo_write(adc1_fifo, signed_buf1, ADC_BUF_SZ/4) == -1) {
            adc_overflow++;
        }
        if (fifo_write(adc2_fifo, signed_buf2, ADC_BUF_SZ/4) == -1) {
            adc_overflow++;
        }

        /* Clear DMA Stream Transfer Complete interrupt pending bit */

        DMA_ClearITPendingBit(DMA2_Stream0, DMA_IT_TCIF0);  
    }
}

