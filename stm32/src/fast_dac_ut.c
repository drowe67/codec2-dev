/*---------------------------------------------------------------------------*\

  FILE........: dac_ut.c
  AUTHOR......: David Rowe
  DATE CREATED: May 31 2013

  Plays a 500 Hz sine wave sampled at 16 kHz out of PA5 on a Discovery board,
  or the speaker output of the SM1000.

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
#include "stm32f4_dacduc.h"

#define SINE_SAMPLES   60


/* 32 sample sine wave which at Fs=16kHz will be 500Hz.  Note samples
   are 16 bit 2's complement, the DAC driver convertsto 12 bit
   unsigned. */

short aWave[] = {4095,0,4095,0,4095,0,4095,0,4095,0,4095,0,4095,0,4095,0,4095,0,4095,0,4095,0,4095,0,4095,0,4095,0,4095,0,4095,0,
	4095,0,4095,0,4095,0,4095,0,4095,0,4095,0,4095,0,4095,0,4095,0,4095,0,4095,0,4095,0,4095,0,4095,0,4095,0,4095,0,};

short aSine[] = {390, 3981, 1410, 841, 4080, 841, 1410, 3981, 390, 2040, 3691, 100, 2671, 3240, 0, 3240, 2671, 100, 3691, 2040, 390, 3981, 1410, 841, 4080, 841, 1410, 3981, 390, 2040, 3691, 100, 2671, 3240, 0, 3240, 2671, 100, 3691, 2040,390, 3981, 1410, 841, 4080, 841, 1410, 3981, 390, 2040, 3691, 100, 2671, 3240, 0, 3240, 2671, 100, 3691, 2040
};

int main(void) {

    dac_open(4*DAC_DUC_BUF_SZ);
    while (1) {

        /* keep DAC FIFOs topped up */

        dac1_write((short*)aWave, SINE_SAMPLES);
    }
   
}

