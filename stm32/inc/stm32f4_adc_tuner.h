/*---------------------------------------------------------------------------*\

  FILE........: stm32f4_adc_tuner.h
  AUTHOR......: David Rowe
  DATE CREATED: 19 Feb 2015

  Single channel ADC driver module for STM32F4 that samples pin PA1 at
  2 MHz and down converts to 50 kHz, with "tuning" centred at 500 kHz.

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

#ifndef __STM32F4_ADC_TUNER__
#define __STM32F4_ADC_TUNER__

#define ADC_TUNER_M  45   /* decimation rate */
#define ADC_TUNER_N  160
#define ADC_TUNER_BUF_SZ  (ADC_TUNER_M*ADC_TUNER_N)

void adc_open(int fifo_sz);
int adc1_read(short buf[], int n); /* ADC1 Pin PA1 */
void adc_set_tuner_en(short flag); /* disable tuner to get raw ADC samples written to fifo */

#endif
