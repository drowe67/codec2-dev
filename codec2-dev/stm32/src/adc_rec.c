/*---------------------------------------------------------------------------*\

  FILE........: adc_rec.c
  AUTHOR......: David Rowe
  DATE CREATED: 30 May 2014

  Recordss a 16 kHz sample rate raw file from the STM32F4 ADC.

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

#include <stdlib.h>
#include "stm32f4_adc.h"
#include "gdb_stdio.h"

#define REC_TIME_SECS 10
#define N   2000
#define FS  16000

int main(void){
    short  buf[N];
    FILE  *frec;
    int    i, bufs;

    adc_open();

    frec = fopen("stm_out.raw", "wb");
    if (frec == NULL) {
        printf("Error opening input file: stm_out.raw\n\nTerminating....\n");
        exit(1);
    }
    bufs = FS*REC_TIME_SECS/N;

    printf("Starting!\n");
    for(i=0; i<bufs; i++) {
        while(adc_read(buf, N) == -1);
        fwrite(buf, sizeof(short), N, frec);  
        printf(".");
    }
    fclose(frec);
    printf("Finished!\n");
}
