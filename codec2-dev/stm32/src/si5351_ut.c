/*---------------------------------------------------------------------------*\

  FILE........: si5351_ut.c
  AUTHOR......: David Rowe
  DATE CREATED: June 2016

  Generates a 10MHz signal on CLK0 ouput of Si5351, should be visible in 
  attenuated form on SP7/SP7 of SM2000.

\*---------------------------------------------------------------------------*/

/*
  Copyright (C) 2016 David Rowe

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
#include "new_i2c.h"
#include "si53xx.h"
//#include "gdb_stdio.h"
//#define printf gdb_stdio_printf

int main(void) {

    //printf("starting\n");
    I2C_Setup();
    si5351_init(0, SI5351_CRYSTAL_LOAD_6PF, 0);
    si5351_set_freq(1000000000ULL, 0, SI5351_CLK0);
    // int ret = si5351_set_freq(1000000000ULL, 0, SI5351_CLK0);
    //printf(" ret = %d\n", ret);

    while(1);
}
