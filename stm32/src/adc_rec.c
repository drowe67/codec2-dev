/*---------------------------------------------------------------------------*\

  FILE........: adc_rec.c
  AUTHOR......: David Rowe
  DATE CREATED: 30 May 2014

  Records a 16 kHz sample rate raw file from one of the ADC channels,
  which are connected to pins PA1 (ADC1) and PA2 (ADC2).

  Note the semi-hosting system isn't fast enough to transfer 2 16 kHz
  streams at once.

  ~/stlink$ sudo ./st-util -f ~/codec2-dev/stm32/adc_rec.elf
  ~/codec2-dev/stm32$ ~/gcc-arm-none-eabi-4_7-2013q1/bin/arm-none-eabi-gdb adc_rec.elf

  (when finished)
  $ play -r 16000 -s -2 ~/stlink/adc.raw

  adc1 -> "from radio"
  adc2 -> "mic amp"

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
#include "stm32f4xx_gpio.h"

#define REC_TIME_SECS 10
#define  N  (ADC_BUF_SZ*6)
#define FS  16000

extern int adc_overflow1;
extern int adc_overflow2;

int main(void){
    short  buf[N];
    FILE  *fadc;
    int    i, bufs;

    fadc = fopen("adc.raw", "wb");
    if (fadc == NULL) {
        printf("Error opening input file: adc.raw\n\nTerminating....\n");
        exit(1);
    }
    bufs = FS*REC_TIME_SECS/N;

    printf("Starting!\n");
    adc_open(ADC_FS_16KHZ, 4*N);

    for(i=0; i<bufs; i++) {
        while(adc2_read(buf, N) == -1);
        fwrite(buf, sizeof(short), N, fadc);
        printf("adc_overflow1: %d  adc_overflow2: %d   \n", adc_overflow1, adc_overflow2);
    }
    fclose(fadc);

    printf("Finished!\n");
}
