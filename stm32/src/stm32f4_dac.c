/*---------------------------------------------------------------------------*\

  FILE........: stm32f4_dac.c
  AUTHOR......: David Rowe
  DATE CREATED: 1 June 2013

  DAC driver module for STM32F4.

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
#include "stm32f4xx.h"
#include "codec2_fifo.h"
#include "stm32f4_dac.h"

#define DAC_DHR12R2_ADDRESS    0x40007414
#define DAC_DHR12L2_ADDRESS    0x40007418

#define DAC_BUF_SZ   320
#define FIFO_SZ      1000
#define DAC_MAX      4096

DAC_InitTypeDef  DAC_InitStructure;
struct FIFO *DMA1_Stream6_fifo;

unsigned short dac_buf[DAC_BUF_SZ];

static void TIM6_Config(void);
static void DAC_Ch2_Config(void);

int dac_underflow;

void dac_open(void) {

    memset(dac_buf, 32768, sizeof(short)*DAC_BUF_SZ);

    /* Create fifo */

    DMA1_Stream6_fifo = fifo_create(FIFO_SZ);
    assert(DMA1_Stream6_fifo != NULL);

    /*!< At this stage the microcontroller clock setting is already configured, 
      this is done through SystemInit() function which is called from startup
      files (startup_stm32f40xx.s/startup_stm32f427x.s) before to branch to 
      application main. 
      To reconfigure the default setting of SystemInit() function, refer to
      system_stm32f4xx.c file
    */    

    /* Preconfiguration before using DAC----------------------------------------*/

    GPIO_InitTypeDef GPIO_InitStructure;

    /* DMA1 clock enable */
    RCC_AHB1PeriphClockCmd(RCC_AHB1Periph_DMA1, ENABLE);
    /* GPIOA clock enable (to be used with DAC) */
    RCC_AHB1PeriphClockCmd(RCC_AHB1Periph_GPIOA, ENABLE);                         
    /* DAC Periph clock enable */
    RCC_APB1PeriphClockCmd(RCC_APB1Periph_DAC, ENABLE);

    /* DAC channel 1 & 2 (DAC_OUT1 = PA.4)(DAC_OUT2 = PA.5) configuration */
    GPIO_InitStructure.GPIO_Pin = GPIO_Pin_4 | GPIO_Pin_5;
    GPIO_InitStructure.GPIO_Mode = GPIO_Mode_AN;
    GPIO_InitStructure.GPIO_PuPd = GPIO_PuPd_NOPULL;
    GPIO_Init(GPIOA, &GPIO_InitStructure);

    /* TIM6 Configuration ------------------------------------------------------*/

    TIM6_Config();  
    DAC_Ch2_Config();
    
}

/* Accepts signed 16 bit samples */

int dac_write(short buf[], int n) {   
    return fifo_write(DMA1_Stream6_fifo, buf, n);
}

/**             
  * @brief  TIM6 Configuration
  * @note   TIM6 configuration is based on APB1 frequency
  * @note   TIM6 Update event occurs each TIM6CLK/256   
  * @param  None
  * @retval None
  */
static void TIM6_Config(void)
{
  TIM_TimeBaseInitTypeDef    TIM_TimeBaseStructure;
  /* TIM6 Periph clock enable */
  RCC_APB1PeriphClockCmd(RCC_APB1Periph_TIM6, ENABLE);
  
  /* --------------------------------------------------------
  
  TIM3 input clock (TIM6CLK) is set to 2 * APB1 clock (PCLK1), since
  APB1 prescaler is different from 1 (see system_stm32f4xx.c and Fig
  13 clock tree figure in DM0031020.pdf).

     Sample rate Fs = 2*PCLK1/TIM_ClockDivision 
                    = (HCLK/2)/TIM_ClockDivision
                    
  ----------------------------------------------------------- */

  /* Time base configuration */
  TIM_TimeBaseStructInit(&TIM_TimeBaseStructure); 
  TIM_TimeBaseStructure.TIM_Period = 5250;          
  TIM_TimeBaseStructure.TIM_Prescaler = 0;       
  TIM_TimeBaseStructure.TIM_ClockDivision = 0;    
  TIM_TimeBaseStructure.TIM_CounterMode = TIM_CounterMode_Up;  
  TIM_TimeBaseInit(TIM6, &TIM_TimeBaseStructure);

  /* TIM6 TRGO selection */

  TIM_SelectOutputTrigger(TIM6, TIM_TRGOSource_Update);
  
  /* TIM6 enable counter */
  TIM_Cmd(TIM6, ENABLE);
}

/**
  * @brief  DAC  Channel2 SineWave Configuration
  * @param  None
  * @retval None
  */
static void DAC_Ch2_Config(void)
{
  DMA_InitTypeDef DMA_InitStructure;
  NVIC_InitTypeDef NVIC_InitStructure;
  
  /* DAC channel2 Configuration */
  DAC_InitStructure.DAC_Trigger = DAC_Trigger_T6_TRGO;
  DAC_InitStructure.DAC_WaveGeneration = DAC_WaveGeneration_None;
  DAC_InitStructure.DAC_OutputBuffer = DAC_OutputBuffer_Disable;
  DAC_Init(DAC_Channel_2, &DAC_InitStructure);

  /* DMA1_Stream6 channel7 configuration **************************************/
  DMA_DeInit(DMA1_Stream6);
  DMA_InitStructure.DMA_Channel = DMA_Channel_7;  
  DMA_InitStructure.DMA_PeripheralBaseAddr = (uint32_t)DAC_DHR12L2_ADDRESS;
  DMA_InitStructure.DMA_Memory0BaseAddr = (uint32_t)dac_buf;
  DMA_InitStructure.DMA_DIR = DMA_DIR_MemoryToPeripheral;
  DMA_InitStructure.DMA_BufferSize = DAC_BUF_SZ;
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
  DMA_Init(DMA1_Stream6, &DMA_InitStructure);

  /* Enable DMA Half & Complete interrupts */
  DMA_ITConfig(DMA1_Stream6, DMA_IT_TC | DMA_IT_HT, ENABLE);

  /* Enable the DMA Stream IRQ Channel */

  NVIC_InitStructure.NVIC_IRQChannel = DMA1_Stream6_IRQn;
  NVIC_InitStructure.NVIC_IRQChannelPreemptionPriority = 0;
  NVIC_InitStructure.NVIC_IRQChannelSubPriority = 0;
  NVIC_InitStructure.NVIC_IRQChannelCmd = ENABLE;
  NVIC_Init(&NVIC_InitStructure);     

  /* Enable DMA1_Stream6 */
  DMA_Cmd(DMA1_Stream6, ENABLE);

  /* Enable DAC Channel2 */
  DAC_Cmd(DAC_Channel_2, ENABLE);

  /* Enable DMA for DAC Channel2 */
  DAC_DMACmd(DAC_Channel_2, ENABLE);
}

/******************************************************************************/
/*                 STM32F4xx Peripherals Interrupt Handlers                   */
/*  Add here the Interrupt Handler for the used peripheral(s) (PPP), for the  */
/*  available peripheral interrupt handler's name please refer to the startup */
/*  file (startup_stm32f40xx.s/startup_stm32f427x.s).                         */
/******************************************************************************/

/*
  This function handles DMA Stream interrupt request.
*/

void DMA1_Stream6_IRQHandler(void) {
    int i, sam;
    short signed_buf[DAC_BUF_SZ/2];

    /* Transfer half empty interrupt */

    if(DMA_GetITStatus(DMA1_Stream6, DMA_IT_HTIF6) != RESET) {
        /* fill first half from fifo */

        if (fifo_read(DMA1_Stream6_fifo, signed_buf, DAC_BUF_SZ/2) == -1) {
            memset(signed_buf, 0, sizeof(short)*DAC_BUF_SZ/2);
            dac_underflow++;
        }

        /* convert to unsigned */

        for(i=0; i<DAC_BUF_SZ/2; i++) {
            sam = (int)signed_buf[i] + 32768;
            dac_buf[i] = (unsigned short)(sam);
        }

        /* Clear DMA Stream Transfer Complete interrupt pending bit */

        DMA_ClearITPendingBit(DMA1_Stream6, DMA_IT_HTIF6);  
    }

    /* Transfer complete interrupt */

    if(DMA_GetITStatus(DMA1_Stream6, DMA_IT_TCIF6) != RESET) {
        /* fill second half from fifo */

        if (fifo_read(DMA1_Stream6_fifo, signed_buf, DAC_BUF_SZ/2) == -1) {
            memset(signed_buf, 0, sizeof(short)*DAC_BUF_SZ/2);
            dac_underflow++;
        }

        /* convert to unsigned */

        for(i=0; i<DAC_BUF_SZ/2; i++) {
            sam = (int)signed_buf[i] + 32768;
            dac_buf[i+DAC_BUF_SZ/2] = (unsigned short)(sam);
        }

        /* Clear DMA Stream Transfer Complete interrupt pending bit */

        DMA_ClearITPendingBit(DMA1_Stream6, DMA_IT_TCIF6);  
    }
}
