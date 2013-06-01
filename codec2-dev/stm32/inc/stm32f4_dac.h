/*---------------------------------------------------------------------------*\

  FILE........: stm32f4_dac.h
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

#ifndef __STM32F4_DAC__
#define __STM32F4_DAC__

void dac_open(void);
int dac_write(short buf[], int n); 

#endif
