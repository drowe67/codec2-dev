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

#define printf gdb_stdio_printf
#define fopen gdb_stdio_fopen
#define fclose gdb_stdio_fclose
#define fread gdb_stdio_fread
#define fwrite gdb_stdio_fwrite

#define N1 24000
#define N2   320

int main(void) {
    short *buf, *pbuf;
    FILE  *fin;
    int    i, nframes;

    buf = (short*)malloc(N1*sizeof(short));
    dac_open();

    fin = fopen("stm_in.raw", "rb");
    if (fin == NULL) {
        printf("Error opening input file: stm_in.raw\n\nTerminating....\n");
        exit(1);
    }
    fread(buf, sizeof(short), N1, fin);
    fclose(fin);

    nframes = N1/N2;
    while(1) {
        printf("Starting!\n");
        pbuf = buf;
        for(i=0; i<nframes; i++) {
            while(dac_write(pbuf, N2) == -1);
            pbuf += N2;
        } 
        printf("Finished!\n");
    }
}

