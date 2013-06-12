/*---------------------------------------------------------------------------*\

  FILE........: stm32f4_adc.c
  AUTHOR......: David Rowe
  DATE CREATED: 4 June 2013

  ADC driver module for STM32F4.

  TODO:
  + just get ADC to run at all, prove its sampling something....
  + timer config to drive ADC conversion, measure sample rate and confirm 16kHz
  + larger ADC DMA buffer
  + half and finished interrupts, ISR
  + fifos
  + work out a way to unit test
 
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
#include "codec2_fifo.h"
#include "gdb_stdio.h"

#include "stm32f4xx_adc.h"
#include "stm32f4xx_gpio.h"
#include "stm32f4xx_rcc.h"
 
#define TRY3
#ifdef TRY3

#define ADCx                     ADC1
#define ADCx_CLK                 RCC_APB2Periph_ADC1
#define DMA_CHANNELx             DMA_Channel_0
#define DMA_STREAMx              DMA2_Stream0
#define ADCx_DR_ADDRESS          ((uint32_t)0x4001204C)

static void ADC_Config(void);

volatile unsigned short uhADCConvertedValue;

int adc_convert(){
 while(!ADC_GetFlagStatus(ADC1, ADC_FLAG_EOC));//Processing the conversion
 return ADC_GetConversionValue(ADC1); //Return the converted data
}

int main(void)
{

  /* ADC1 Channel Vbat configuration */
  ADC_Config();
 
  /* Start ADC1 Software Conversion */ 
  ADC_SoftwareStartConv(ADC1);

  while (1) {
      //uhADCConvertedValue = adc_convert();
      printf("try 3: %d\n", uhADCConvertedValue);
  }
}

static void ADC_Config(void);

/**
  * @brief  ADC1 Channel Vbat configuration
  * @note   This function Configure the ADC peripheral  
            1) Enable peripheral clocks
            2) DMA2_Stream0 channel 0 configuration
            3) Configure ADC1 Channel18 (VBAT) 
  * @param  None
  * @retval None
  */

static void ADC_Config(void)
{
  ADC_InitTypeDef       ADC_InitStructure;
  ADC_CommonInitTypeDef ADC_CommonInitStructure;
  DMA_InitTypeDef       DMA_InitStructure;

    
  /* Enable peripheral clocks *************************************************/
  RCC_AHB1PeriphClockCmd(RCC_AHB1Periph_DMA2, ENABLE);
  RCC_APB2PeriphClockCmd(ADCx_CLK, ENABLE);

  /* DMA2_Stream0 channel0 configuration **************************************/
  DMA_DeInit(DMA2_Stream0);
  DMA_InitStructure.DMA_Channel = DMA_CHANNELx;  
  DMA_InitStructure.DMA_PeripheralBaseAddr = (uint32_t)ADCx_DR_ADDRESS;
  DMA_InitStructure.DMA_Memory0BaseAddr = (uint32_t)&uhADCConvertedValue;
  DMA_InitStructure.DMA_DIR = DMA_DIR_PeripheralToMemory;
  DMA_InitStructure.DMA_BufferSize = 1;
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

#define TMP
#ifdef TMP
  /* ADC Common Init **********************************************************/
  ADC_CommonInitStructure.ADC_Mode = ADC_Mode_Independent;
  ADC_CommonInitStructure.ADC_Prescaler = ADC_Prescaler_Div2;
  ADC_CommonInitStructure.ADC_DMAAccessMode = ADC_DMAAccessMode_Disabled;
  ADC_CommonInitStructure.ADC_TwoSamplingDelay = ADC_TwoSamplingDelay_5Cycles;
  ADC_CommonInit(&ADC_CommonInitStructure);
#endif

  /* ADC1 Init ****************************************************************/
  ADC_InitStructure.ADC_Resolution = ADC_Resolution_12b;
  ADC_InitStructure.ADC_ScanConvMode = DISABLE;
  ADC_InitStructure.ADC_ContinuousConvMode = ENABLE;
  ADC_InitStructure.ADC_ExternalTrigConvEdge = ADC_ExternalTrigConvEdge_None;
  ADC_InitStructure.ADC_ExternalTrigConv = ADC_ExternalTrigConv_T1_CC1;
  ADC_InitStructure.ADC_DataAlign = ADC_DataAlign_Right;
  ADC_InitStructure.ADC_NbrOfConversion = 1;
  ADC_Init(ADCx, &ADC_InitStructure);

  //#define VBAT
#ifdef VBAT
  /* ADC1 regular channel18 (VBAT) configuration ******************************/
  ADC_RegularChannelConfig(ADCx, ADC_Channel_Vbat, 1, ADC_SampleTime_15Cycles);

  /* Enable VBAT channel */
  ADC_VBATCmd(ENABLE); 
#else
  ADC_RegularChannelConfig(ADC1, ADC_Channel_TempSensor, 1, ADC_SampleTime_15Cycles);
  ADC_TempSensorVrefintCmd(ENABLE);
#endif

  /* Enable DMA request after last transfer (Single-ADC mode) */
  ADC_DMARequestAfterLastTransferCmd(ADCx, ENABLE);

  /* Enable ADC1 DMA */
  ADC_DMACmd(ADCx, ENABLE);
  
  /* DMA2_Stream0 enable */
  DMA_Cmd(DMA_STREAMx, ENABLE);
    

  /* Enable ADC1 **************************************************************/
  ADC_Cmd(ADCx, ENABLE);
}
#endif

//#define TRY2
#ifdef TRY2

int ConvertedValue = 0; //Converted value readed from ADC
#define ADCx_DR_ADDRESS          ((uint32_t)0x4001204C)
volatile unsigned short uhADCConvertedValue;

#define DMA_CHANNELx             DMA_Channel_0
#define DMA_STREAMx              DMA2_Stream0
#define ADCx                     ADC1


void adc_configure(){
    ADC_InitTypeDef ADC_init_structure; //Structure for adc confguration
    GPIO_InitTypeDef GPIO_initStructre; //Structure for analog input pin
    DMA_InitTypeDef       DMA_InitStructure;

    // Clock configuration

    RCC_APB2PeriphClockCmd(RCC_APB2Periph_ADC1,ENABLE);
    RCC_AHB1PeriphClockCmd(RCC_AHB1ENR_GPIOCEN,ENABLE);
    RCC_AHB1PeriphClockCmd(RCC_AHB1Periph_DMA2, ENABLE);

    // Analog pin configuration

    GPIO_initStructre.GPIO_Pin = GPIO_Pin_0;//The channel 10 is connected to PC0
    GPIO_initStructre.GPIO_Mode = GPIO_Mode_AN; //The PC0 pin is configured in analog mode
    GPIO_initStructre.GPIO_PuPd = GPIO_PuPd_NOPULL; //We don't need any pull up or pull down
    GPIO_Init(GPIOC,&GPIO_initStructre);//Affecting the port with the initialization structure configuration

    // ADC structure configuration

    ADC_DeInit();
    ADC_init_structure.ADC_DataAlign = ADC_DataAlign_Right;//data converted will be shifted to right
    ADC_init_structure.ADC_Resolution = ADC_Resolution_12b;//Input voltage is converted into a 12bit number giving a maximum value of 4096
    ADC_init_structure.ADC_ContinuousConvMode = ENABLE; //the conversion is continuous, the input data is converted more than once
    ADC_init_structure.ADC_ExternalTrigConv = ADC_ExternalTrigConv_T1_CC1;// conversion is synchronous with TIM1 and CC1 (actually I'm not sure about this one :/)
    ADC_init_structure.ADC_ExternalTrigConvEdge = ADC_ExternalTrigConvEdge_None;//no trigger for conversion
    ADC_init_structure.ADC_NbrOfConversion = 1;//I think this one is clear :p
    ADC_init_structure.ADC_ScanConvMode = DISABLE;//The scan is configured in one channel
    ADC_Init(ADCx,&ADC_init_structure);//Initialize ADC with the previous configuration

    // Select the channel to be read from

    ADC_RegularChannelConfig(ADCx,ADC_Channel_10,1,ADC_SampleTime_144Cycles);

    /* DMA  configuration **************************************/

    DMA_DeInit(DMA_STREAMx);
    DMA_InitStructure.DMA_Channel = DMA_CHANNELx;  
    DMA_InitStructure.DMA_PeripheralBaseAddr = (uint32_t)ADCx_DR_ADDRESS;
    DMA_InitStructure.DMA_Memory0BaseAddr = (uint32_t)&uhADCConvertedValue;
    DMA_InitStructure.DMA_DIR = DMA_DIR_PeripheralToMemory;
    DMA_InitStructure.DMA_BufferSize = 1;
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

    // Enable ADC conversion

    ADC_Cmd(ADC1,ENABLE);

}

int adc_convert(){
    while(!ADC_GetFlagStatus(ADC1, ADC_FLAG_EOC));//Processing the conversion
    return ADC_GetConversionValue(ADC1); //Return the converted data
}

int main(void){
    adc_configure();
    ADC_SoftwareStartConv(ADC1);
    while(1){
        //ConvertedValue = adc_convert();
        //printf("ConvertedValue = %d\n", ConvertedValue); 
        printf("Try 2: %d\n",uhADCConvertedValue); 
    }
}
#endif

#ifdef TRY1

#define ADCx                     ADC1
#define ADC_CHANNEL              ADC_Channel_7
#define ADCx_CLK                 RCC_APB2Periph_ADC1
#define ADCx_CHANNEL_GPIO_CLK    RCC_AHB1Periph_GPIOA
#define GPIO_PIN                 GPIO_Pin_7
#define GPIO_PORT                GPIOA
#define DMA_CHANNELx             DMA_Channel_0
#define DMA_STREAMx              DMA2_Stream0
#define ADCx_DR_ADDRESS          ((uint32_t)0x4001204C)

#define TIM1_CCR3_ADDRESS    0x4001223C

TIM_TimeBaseInitTypeDef  TIM_TimeBaseStructure;
TIM_OCInitTypeDef  TIM_OCInitStructure;

uint16_t uhTimerPeriod;
uint16_t aSRC_Buffer[3] = {0, 0, 0};

unsigned short uhADCxConvertedValue;

static void ADC_Config(void);
void Timer1Config();
static void TIM_Config(void);

int main(void) {
    /* check ADC value changes with GND and VCC sampled */
    /* look for 17kHz (ish) on timer pins */

    ADC_Config();
    ADC_SoftwareStartConv(ADCx);

    Timer1Config();

    /* todo: that start up cmd here */

    while(1) {
        printf("ADC: %d\n", uhADCxConvertedValue);
    }

}

/**
  * @brief  ADC3 channel07 with DMA configuration
  * @note   This function Configure the ADC peripheral  
            1) Enable peripheral clocks
            2) DMA2_Stream0 channel2 configuration
            3) Configure ADC Channel7 pin as analog input
            4) Configure ADC3 Channel7 
  * @param  None
  * @retval None
  */
static void ADC_Config(void)
{
  ADC_InitTypeDef       ADC_InitStructure;
  ADC_CommonInitTypeDef ADC_CommonInitStructure;
  DMA_InitTypeDef       DMA_InitStructure;
  GPIO_InitTypeDef      GPIO_InitStructure;

  /* Enable ADCx, DMA and GPIO clocks ****************************************/ 
  RCC_AHB1PeriphClockCmd(RCC_AHB1Periph_DMA2, ENABLE);
  RCC_AHB1PeriphClockCmd(ADCx_CHANNEL_GPIO_CLK, ENABLE);  
  RCC_APB2PeriphClockCmd(ADCx_CLK, ENABLE);
  

  /* DMA2 Stream0 channel2 configuration **************************************/
  DMA_InitStructure.DMA_Channel = DMA_CHANNELx;  
  DMA_InitStructure.DMA_PeripheralBaseAddr = (uint32_t)ADCx_DR_ADDRESS;
  DMA_InitStructure.DMA_Memory0BaseAddr = (uint32_t)&uhADCxConvertedValue;
  DMA_InitStructure.DMA_DIR = DMA_DIR_PeripheralToMemory;
  DMA_InitStructure.DMA_BufferSize = 1;
  DMA_InitStructure.DMA_PeripheralInc = DMA_PeripheralInc_Disable;
  DMA_InitStructure.DMA_MemoryInc = DMA_MemoryInc_Disable;
  DMA_InitStructure.DMA_PeripheralDataSize = DMA_PeripheralDataSize_HalfWord;
  DMA_InitStructure.DMA_MemoryDataSize = DMA_MemoryDataSize_HalfWord;
  DMA_InitStructure.DMA_Mode = DMA_Mode_Circular;
  DMA_InitStructure.DMA_Priority = DMA_Priority_High;
  DMA_InitStructure.DMA_FIFOMode = DMA_FIFOMode_Disable;         
  DMA_InitStructure.DMA_FIFOThreshold = DMA_FIFOThreshold_HalfFull;
  DMA_InitStructure.DMA_MemoryBurst = DMA_MemoryBurst_Single;
  DMA_InitStructure.DMA_PeripheralBurst = DMA_PeripheralBurst_Single;
  DMA_Init(DMA_STREAMx, &DMA_InitStructure);
  DMA_Cmd(DMA_STREAMx, ENABLE);

  /* Configure ADC3 Channel7 pin as analog input ******************************/
  GPIO_InitStructure.GPIO_Pin = GPIO_PIN;
  GPIO_InitStructure.GPIO_Mode = GPIO_Mode_AN;
  GPIO_InitStructure.GPIO_PuPd = GPIO_PuPd_NOPULL ;
  GPIO_Init(GPIO_PORT, &GPIO_InitStructure);

  /* ADC Common Init **********************************************************/
  ADC_CommonInitStructure.ADC_Mode = ADC_Mode_Independent;
  ADC_CommonInitStructure.ADC_Prescaler = ADC_Prescaler_Div2;
  ADC_CommonInitStructure.ADC_DMAAccessMode = ADC_DMAAccessMode_Disabled;
  ADC_CommonInitStructure.ADC_TwoSamplingDelay = ADC_TwoSamplingDelay_5Cycles;
  ADC_CommonInit(&ADC_CommonInitStructure);

  /* ADC3 Init ****************************************************************/
  ADC_InitStructure.ADC_Resolution = ADC_Resolution_12b;
  ADC_InitStructure.ADC_ScanConvMode = DISABLE;
  ADC_InitStructure.ADC_ContinuousConvMode = ENABLE;
  ADC_InitStructure.ADC_ExternalTrigConvEdge = ADC_ExternalTrigConvEdge_None;
  ADC_InitStructure.ADC_ExternalTrigConv = ADC_ExternalTrigConv_T1_CC1;
  ADC_InitStructure.ADC_DataAlign = ADC_DataAlign_Right;
  ADC_InitStructure.ADC_NbrOfConversion = 1;
  ADC_Init(ADCx, &ADC_InitStructure);

  /* ADC3 regular channel7 configuration *************************************/
  ADC_RegularChannelConfig(ADCx, ADC_CHANNEL, 1, ADC_SampleTime_3Cycles);

  /* Enable DMA request after last transfer (Single-ADC mode) */
  ADC_DMARequestAfterLastTransferCmd(ADCx, ENABLE);

  /* Enable ADC3 DMA */
  ADC_DMACmd(ADCx, ENABLE);

  /* Enable ADC3 */
  ADC_Cmd(ADCx, ENABLE);
}

void Timer1Config() {

    /* TIM Configuration */

    TIM_Config();

    /* TIM1 DMA Transfer example -------------------------------------------------
  
       TIM1 input clock (TIM1CLK) is set to 2 * APB2 clock (PCLK2), since APB2 
       prescaler is different from 1.   
       TIM1CLK = 2 * PCLK2  
       PCLK2 = HCLK / 2 
       => TIM1CLK = 2 * (HCLK / 2) = HCLK = SystemCoreClock
  
       TIM1CLK = SystemCoreClock, Prescaler = 0, TIM1 counter clock = SystemCoreClock
       SystemCoreClock is set to 168 MHz for STM32F4xx devices.

       The objective is to configure TIM1 channel 3 to generate complementary PWM
       signal with a frequency equal to 17.57 KHz:
       - TIM1_Period = (SystemCoreClock / 17570) - 1
       and a variable duty cycle that is changed by the DMA after a specific number of
       Update DMA request.

       The number of this repetitive requests is defined by the TIM1 Repetion counter,
       each 3 Update Requests, the TIM1 Channel 3 Duty Cycle changes to the next new 
       value defined by the aSRC_Buffer.
  
       Note: 
       SystemCoreClock variable holds HCLK frequency and is defined in system_stm32f4xx.c file.
       Each time the core clock (HCLK) changes, user had to call SystemCoreClockUpdate()
       function to update SystemCoreClock variable value. Otherwise, any configuration
       based on this variable will be incorrect.  
       -----------------------------------------------------------------------------*/
  
    /* Compute the value to be set in ARR regiter to generate signal frequency at 17.57 Khz */
    uhTimerPeriod = (SystemCoreClock / 17570 ) - 1;
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
    TIM_TimeBaseStructure.TIM_Prescaler = 0;
    TIM_TimeBaseStructure.TIM_CounterMode = TIM_CounterMode_Up;
    TIM_TimeBaseStructure.TIM_Period = uhTimerPeriod;
    TIM_TimeBaseStructure.TIM_ClockDivision = 0;
    TIM_TimeBaseStructure.TIM_RepetitionCounter = 3;

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

/* DR: note I dont think any of this is needed, except perhaps to check
   timer frequency */

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
  DMA_InitStructure.DMA_BufferSize = 3;
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
#endif
