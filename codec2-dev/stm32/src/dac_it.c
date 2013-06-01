/**
  ******************************************************************************
  * @file    DMA/DMA_FLASHToRAM/stm32f4xx_it.c 
  * @author  MCD Application Team
  * @version V1.1.0
  * @date    18-January-2013
  * @brief   Main Interrupt Service Routines.
  *          This file provides template for all exceptions handler and 
  *          peripherals interrupt service routine.
  ******************************************************************************
  * @attention
  *
  * <h2><center>&copy; COPYRIGHT 2013 STMicroelectronics</center></h2>
  *
  * Licensed under MCD-ST Liberty SW License Agreement V2, (the "License");
  * You may not use this file except in compliance with the License.
  * You may obtain a copy of the License at:
  *
  *        http://www.st.com/software_license_agreement_liberty_v2
  *
  * Unless required by applicable law or agreed to in writing, software 
  * distributed under the License is distributed on an "AS IS" BASIS, 
  * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
  * See the License for the specific language governing permissions and
  * limitations under the License.
  *
  ******************************************************************************
  */
int interrupts;


/* Includes ------------------------------------------------------------------*/
#include "dac_it.h"

/** @addtogroup STM32F4xx_StdPeriph_Examples
  * @{
  */

/** @addtogroup DMA_FLASHToRAM
  * @{
  */  

/* Private typedef -----------------------------------------------------------*/
/* Private define ------------------------------------------------------------*/
/* Private macro -------------------------------------------------------------*/
/* Private variables ---------------------------------------------------------*/
/* Private function prototypes -----------------------------------------------*/
/* Private functions ---------------------------------------------------------*/

/******************************************************************************/
/*            Cortex-M4 Processor Exceptions Handlers                         */
/******************************************************************************/

/**
  * @brief   This function handles NMI exception.
  * @param  None
  * @retval None
  */
void NMI_Handler(void)
{
}

/**
  * @brief  This function handles Hard Fault exception.
  * @param  None
  * @retval None
  */
void HardFault_Handler(void)
{
  /* Go to infinite loop when Hard Fault exception occurs */
  while (1)
  {
  }
}

/**
  * @brief  This function handles Memory Manage exception.
  * @param  None
  * @retval None
  */
void MemManage_Handler(void)
{
  /* Go to infinite loop when Memory Manage exception occurs */
  while (1)
  {
  }
}

/**
  * @brief  This function handles Bus Fault exception.
  * @param  None
  * @retval None
  */
void BusFault_Handler(void)
{
  /* Go to infinite loop when Bus Fault exception occurs */
  while (1)
  {
  }
}

/**
  * @brief  This function handles Usage Fault exception.
  * @param  None
  * @retval None
  */
void UsageFault_Handler(void)
{
  /* Go to infinite loop when Usage Fault exception occurs */
  while (1)
  {
  }
}

/**
  * @brief  This function handles SVCall exception.
  * @param  None
  * @retval None
  */
void SVC_Handler(void)
{
}

/**
  * @brief  This function handles Debug Monitor exception.
  * @param  None
  * @retval None
  */
void DebugMon_Handler(void)
{
}

/**
  * @brief  This function handles PendSVC exception.
  * @param  None
  * @retval None
  */
void PendSV_Handler(void)
{
}

/**
  * @brief  This function handles SysTick Handler.
  * @param  None
  * @retval None
  */
void SysTick_Handler(void)
{
}

/******************************************************************************/
/*                 STM32F4xx Peripherals Interrupt Handlers                   */
/*  Add here the Interrupt Handler for the used peripheral(s) (PPP), for the  */
/*  available peripheral interrupt handler's name please refer to the startup */
/*  file (startup_stm32f40xx.s/startup_stm32f427x.s).                         */
/******************************************************************************/

/**
  * @brief  This function handles DMA Stream interrupt request.
  * @param  None
  * @retval None
  */
void DMA1_Stream6_IRQHandler(void)
{

  /* Transfer half empty interrupt */

  if(DMA_GetITStatus(DMA1_Stream6, DMA_IT_HTIF6) != RESET))
  {
      /* fill first half from fifo */

      fifo_read(DMA1_Stream6_fifo, dac_buf, DAC_BUF_SZ/2);

      /* Clear DMA Stream Transfer Complete interrupt pending bit */

      DMA_ClearITPendingBit(DMA1_Stream6, DMA_IT_HTIF6);  

      interrupts++;
  }

  /* Transfer complete interrupt */

  if(DMA_GetITStatus(DMA1_Stream6, DMA_IT_TCIF6) != RESET))
  {
      /* fill second half from fifo */

      fifo_read(DMA1_Stream6_fifo, &dac_buf[DAC_BUF_SZ/2], DAC_BUF_SZ/2);

      /* Clear DMA Stream Transfer Complete interrupt pending bit */

      DMA_ClearITPendingBit(DMA1_Stream6, DMA_IT_TCIF6);  

      interrupts++;
  }
}

/**
  * @}
  */ 

/**
  * @}
  */

/************************ (C) COPYRIGHT STMicroelectronics *****END OF FILE****/
