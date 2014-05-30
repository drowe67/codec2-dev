/*---------------------------------------------------------------------------*\

  FILE........: stm32f4_adc.h
  AUTHOR......: David Rowe
  DATE CREATED: 30 May 2014

  ADC driver module for STM32F4.

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

#ifndef __STM32F4_ADC__
#define __STM32F4_ADC__

void adc_open(void);
int adc_read(short buf[], int n); 

#endif
