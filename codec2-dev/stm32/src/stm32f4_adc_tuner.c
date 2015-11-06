/*---------------------------------------------------------------------------*\

  FILE........: stm32f4_adc_tuner.c
  AUTHOR......: David Rowe
  DATE CREATED: 19 Feb 2015

  Single channel ADC driver module for STM32F4 that samples pin PA1 at
  2 MHz and down converts to 50 kHz, with "tuning" centred at 500 kHz.

  See codec2-dev/octave.m for a simulation model.

\*---------------------------------------------------------------------------*/

/*
  Copyright (C) 2015 David Rowe

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
#include "stm32f4_adc_tuner.h"
#include "debugblinky.h"
#include "iir_tuner.h"

struct FIFO *adc1_fifo;
unsigned short adc_buf[ADC_TUNER_BUF_SZ], *padc_buf;
int adc_overflow1;
int half,full;
static short tuner_en = 1;

#define ADCx_DR_ADDRESS          ((uint32_t)0x4001204C)
#define DMA_CHANNELx             DMA_Channel_0
#define DMA_STREAMx              DMA2_Stream0
#define ADCx                     ADC1

void adc_configure();

static void tim2_config(void);

//#define DUMMY_SIGNAL
#ifdef DUMMY_SIGNAL
unsigned short sine[ADC_TUNER_BUF_SZ];
#endif

void adc_open(int fifo_sz) {
    adc1_fifo = fifo_create(fifo_sz);
    assert(adc1_fifo != NULL);

    tim2_config();
    adc_configure();
    init_debug_blinky();
}


/* n signed 16 bit samples in buf[] if return != -1 */

int adc1_read(short buf[], int n) {
    return fifo_read(adc1_fifo, buf, n);
}


void adc_set_tuner_en(short flag)
{
    tuner_en = flag;
}

static void tim2_config(void)
{
  TIM_TimeBaseInitTypeDef    TIM_TimeBaseStructure;

  /* TIM2 Periph clock enable */
  RCC_APB1PeriphClockCmd(RCC_APB1Periph_TIM2, ENABLE);

  /* --------------------------------------------------------

  TIM2 input clock (TIM2CLK) is set to 2 * APB1 clock (PCLK1), since
  APB1 prescaler is different from 1 (see system_stm32f4xx.c and Fig
  13 clock tree figure in DM0031020.pdf).

     Sample rate Fs = 2*PCLK1/)TIM_ClockDivision+1)
                    = (HCLK/2)/(TIM_ClockDivision+1)

  Note from David: The +1 was discovered empirically, still not sure
  if it's right.

  ----------------------------------------------------------- */

  /* Time base configuration */

  TIM_TimeBaseStructInit(&TIM_TimeBaseStructure);
  TIM_TimeBaseStructure.TIM_Period = 41;
  TIM_TimeBaseStructure.TIM_Prescaler = 0;
  TIM_TimeBaseStructure.TIM_ClockDivision = 0;
  TIM_TimeBaseStructure.TIM_CounterMode = TIM_CounterMode_Up;
  TIM_TimeBaseInit(TIM2, &TIM_TimeBaseStructure);

  /* TIM2 TRGO selection */

  TIM_SelectOutputTrigger(TIM2, TIM_TRGOSource_Update);

  /* TIM2 enable counter */

  TIM_Cmd(TIM2, ENABLE);
}


void adc_configure() {
    ADC_InitTypeDef  ADC_init_structure;
    GPIO_InitTypeDef GPIO_initStructre;
    DMA_InitTypeDef  DMA_InitStructure;
    NVIC_InitTypeDef NVIC_InitStructure;

    // Clock configuration

    RCC_APB2PeriphClockCmd(RCC_APB2Periph_ADC1,ENABLE);
    RCC_AHB1PeriphClockCmd(RCC_AHB1ENR_GPIOAEN,ENABLE);
    RCC_AHB1PeriphClockCmd(RCC_AHB1Periph_DMA2, ENABLE);

    // Analog pin configuration ADC1->PA1

    GPIO_initStructre.GPIO_Pin =  GPIO_Pin_1;
    GPIO_initStructre.GPIO_Mode = GPIO_Mode_AN;
    GPIO_initStructre.GPIO_PuPd = GPIO_PuPd_NOPULL;
    GPIO_Init(GPIOA,&GPIO_initStructre);

    // ADC structure configuration

    ADC_DeInit();
    ADC_init_structure.ADC_DataAlign = ADC_DataAlign_Left;
    ADC_init_structure.ADC_Resolution = ADC_Resolution_12b;
    ADC_init_structure.ADC_ContinuousConvMode = DISABLE;
    ADC_init_structure.ADC_ExternalTrigConv = ADC_ExternalTrigConv_T2_TRGO;
    ADC_init_structure.ADC_ExternalTrigConvEdge = ADC_ExternalTrigConvEdge_Rising;
    ADC_init_structure.ADC_NbrOfConversion = 1;
    ADC_Init(ADCx,&ADC_init_structure);

    // Select the channel to be read from

    ADC_RegularChannelConfig(ADCx,ADC_Channel_1,1,ADC_SampleTime_3Cycles);

    /* DMA  configuration **************************************/

    DMA_DeInit(DMA_STREAMx);
    DMA_InitStructure.DMA_Channel = DMA_CHANNELx;
    DMA_InitStructure.DMA_PeripheralBaseAddr = (uint32_t)ADCx_DR_ADDRESS;
    DMA_InitStructure.DMA_Memory0BaseAddr = (uint32_t)adc_buf;
    DMA_InitStructure.DMA_DIR = DMA_DIR_PeripheralToMemory;
    DMA_InitStructure.DMA_BufferSize = ADC_TUNER_BUF_SZ;
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

    // Enable and start ADC conversion

    ADC_Cmd(ADC1,ENABLE);
    ADC_SoftwareStartConv(ADC1);

    padc_buf = adc_buf;

    #ifdef DUMMY_SIGNAL
    int i;

    /* Fs/4 sine wave, right in the middle of the pass band ! */

    for(i=0; i<ADC_TUNER_BUF_SZ; i++)
        sine[i] = 32767;
    for(i=1; i<ADC_TUNER_BUF_SZ; i+=4)
        sine[i] += 32767/4;
    for(i=3; i<ADC_TUNER_BUF_SZ; i+=4)
        sine[i] -= 32767/4;
    padc_buf = sine;
    #endif

}


/*
  This function handles DMA Stream interrupt request.

  ADC_TUNER_BUF_SZ = 45 * 160 = 7200, so one interrupt every 7200/2 = 3600 samples
  or interrupts at a rate of 2E6/3600 = 555.56 Hz.
*/

void DMA2_Stream0_IRQHandler(void) {
    float dec_buf[ADC_TUNER_N/2];

    /* PE0 is asserted high for the duration of this ISR */

    GPIOE->ODR |= (1 << 0);

    /* Half transfer interrupt */

    if(DMA_GetITStatus(DMA2_Stream0, DMA_IT_HTIF0) != RESET) {
        half++;

        if (tuner_en) {
            iir_tuner(dec_buf, padc_buf);

            /* write first half to fifo.  Note we are writing ADC_TUNER_N/2 floats,
               which is equivalent to ADC_TUNER_N shorts.  */

           if (fifo_write(adc1_fifo, (short*)dec_buf, ADC_TUNER_N) == -1) {
                adc_overflow1++;
            }
        }
        else // note: we dump signed shorts when tuner off
            fifo_write(adc1_fifo, (short*)padc_buf, ADC_TUNER_BUF_SZ/2); 

        /* Clear DMA Stream Transfer Complete interrupt pending bit */

        DMA_ClearITPendingBit(DMA2_Stream0, DMA_IT_HTIF0);
    }

    /* Transfer complete interrupt */

    if(DMA_GetITStatus(DMA2_Stream0, DMA_IT_TCIF0) != RESET) {
        full++;

        if (tuner_en) {
            iir_tuner(dec_buf, &padc_buf[ADC_TUNER_BUF_SZ/2]);

            /* write second half to fifo */

            if (fifo_write(adc1_fifo, (short*)dec_buf, ADC_TUNER_N) == -1) {
              adc_overflow1++;
            }
        }
        else
            fifo_write(adc1_fifo, (short*)&padc_buf[ADC_TUNER_BUF_SZ/2], ADC_TUNER_BUF_SZ/2);

        /* Clear DMA Stream Transfer Complete interrupt pending bit */

        DMA_ClearITPendingBit(DMA2_Stream0, DMA_IT_TCIF0);
    }

    GPIOE->ODR &= ~(1 << 0);
}

