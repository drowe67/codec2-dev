/*---------------------------------------------------------------------------*\

  FILE........: iir_tuner.h
  AUTHOR......: David Rowe
  DATE CREATED: 20 Feb 2015

  Header file for IIR tuner function.

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

#ifndef __IIR_TUNER__
#define __IIR_TUNER__

#define IIR_TUNER_DEC_50_10_FILT_MEM 100

void iir_tuner(float dec_50[], unsigned short adc_buf[]);
void iir_tuner_dec_50_to_10(float dec_10[], float dec_50[], int n);

#endif
