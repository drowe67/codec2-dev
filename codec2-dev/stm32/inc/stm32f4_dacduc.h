/*---------------------------------------------------------------------------*\

  FILE........: stm32f4_dac.h
  AUTHOR......: David Rowe
  DATE CREATED: 1 June 2013

  Two channel FIFO buffered DAC driver module for STM32F4. DAC1 is fixed at
  Fs=2Mhz

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

#ifndef __STM32F4_DAC__
#define __STM32F4_DAC__

#define DUC_N 160
#define DUC_M  25
#define DUC_48N 96                     //This is 3/5th DUC_N
#define DAC_DUC_BUF_SZ DUC_M*DUC_N
#define DAC_BUF_SZ 2048

void fast_dac_open(int dac1_fifo_size,int dac2_fifo_size);
int dac1_write(short buf[], int n); /* DAC1 pin PA4 */
int dac2_write(short buf[], int n); /* DAC2 pin PA5 */

#endif
