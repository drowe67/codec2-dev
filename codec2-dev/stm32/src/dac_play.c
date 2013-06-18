/*---------------------------------------------------------------------------*\

  FILE........: dac_play.c
  AUTHOR......: David Rowe
  DATE CREATED: 1 June 2013

  Plays a 16 kHz sample rate raw file to the Discovery DAC.

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

#include <stdlib.h>
#include "stm32f4_dac.h"
#include "gdb_stdio.h"

#define N   2000

int main(void) {
    short  buf[N];
    FILE  *fplay;

    dac_open();

    while(1) {
        fplay = fopen("stm_in.raw", "rb");
        if (fplay == NULL) {
            printf("Error opening input file: stm_in.raw\n\nTerminating....\n");
            exit(1);
        }
    
        printf("Starting!\n");

        while(fread(buf, sizeof(short), N, fplay) == N) {
            while(dac_write(buf, N) == -1);
        }  

        printf("Finished!\n");
        fclose(fplay);
    }

    /* let FIFO empty */

    while(1);
}

