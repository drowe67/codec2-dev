/*---------------------------------------------------------------------------*\

  FILE........: adc_sd.c
  AUTHOR......: David Rowe
  DATE CREATED: 30 May 2014

  Measures the std deviation of the ADC signals.  Used to check noise
  levels on each ADC.

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
#include <math.h>
#include "stm32f4_adc.h"
#include "stm32f4_dac.h"
#include "gdb_stdio.h"

#define REC_TIME_SECS 10
#define  N  (ADC_BUF_SZ*4)
#define FS  16000

static float calc_sd(short x[], int n) {
    float sum, mean, sum_diff, sd;
    int   i;

    sum = 0.0;
    for(i=0; i<n;i++) {
        sum += (float)x[i];
    }
    mean = sum/n;

    sum_diff = 0.0;
    for(i=0; i<n;i++) {
        sum_diff += ((float)x[i] - mean)*((float)x[i] - mean);
    }

    sd = sqrtf(sum_diff/n);

    return sd;
}

int main(void){
    short  buf[N];
    float  sd1, sd2;

    adc_open(ADC_FS_16KHZ, 2*N);

    printf("Starting!\n");
    while(1) {
        while(adc1_read(buf, N) == -1);
        sd1 = calc_sd(buf, N);
        while(adc2_read(buf, N) == -1);
        sd2 = calc_sd(buf, N);

        printf("adc1: %5.1f adc2: %5.1f\n", (double)sd1, (double)sd2);
    }

}
