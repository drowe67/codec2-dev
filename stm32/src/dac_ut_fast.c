/*---------------------------------------------------------------------------*\

  FILE........: dac_ut_fast.c
  AUTHOR......: David Rowe
  DATE CREATED: Sep 2015

  Plays a Fs/4 sine wave sampled out of PA5 on a Discovery board, used for 
  testing high speed DAC operation, e.g. for IF/RF generation.

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
#include "stm32f4_dac.h"


int main(void) {
    dac_open(4*DAC_BUF_SZ);
    while (1);
}

