/*---------------------------------------------------------------------------*\

  FILE........: power_ut.c
  AUTHOR......: David Rowe
  DATE CREATED: 30 May 2014

  Runs Codec 2, ADC, and DAC, to fully exercise STM32C so we can a feel for
  run-time power consumption for SM1000 and hence dimension regulators.

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

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>

#include "stm32f4xx_conf.h"
#include "stm32f4xx.h"
#include "stm32f4_adc.h"
#include "stm32f4_dac.h"
#include "gdb_stdio.h"
#include "codec2.h"
#include "dump.h"
#include "sine.h"
#include "machdep.h"

#ifdef __EMBEDDED__
#define printf gdb_stdio_printf
#define fopen gdb_stdio_fopen
#define fclose gdb_stdio_fclose
#define fread gdb_stdio_fread
#define fwrite gdb_stdio_fwrite
#endif

#define SPEED_TEST_SAMPLES 24000

/* modification of test used to measure codec2 execuation speed.  We read/write ADC/DAC
   but dont do anything with the samples, as they are at 16 kHz and codec needs 8 kHz.  Just
   trying to exercise everything to get a feel for power consumption */

static void c2speedtest(int mode, char inputfile[])
{
    struct CODEC2 *codec2;
    short         *inbuf, *outbuf, *pinbuf, *dummy_buf;
    unsigned char *bits;
    int            nsam, nbit, nframes;
    FILE          *fin;
    int            f, nread;

    codec2 = codec2_create(mode);
    nsam = codec2_samples_per_frame(codec2);
    nframes = SPEED_TEST_SAMPLES/nsam;
    outbuf = (short*)malloc(nsam*sizeof(short));
    inbuf = (short*)malloc(SPEED_TEST_SAMPLES*sizeof(short));
    nbit = codec2_bits_per_frame(codec2);
    bits = (unsigned char*)malloc(nbit*sizeof(char));
    dummy_buf = (short*)malloc(2*nsam*sizeof(short));

    fin = fopen(inputfile, "rb");
    if (fin == NULL) {
        printf("Error opening input file: %s\nTerminating....\n",inputfile);
        exit(1);
    }

    printf("reading samples ....\n");
    nread = fread(inbuf, sizeof(short), SPEED_TEST_SAMPLES, fin);
    if (nread != SPEED_TEST_SAMPLES) {
        printf("error reading %s, %d samples reqd, %d read\n",
               inputfile, SPEED_TEST_SAMPLES, nread);
    }
    fclose(fin);

    pinbuf = inbuf;
    for(f=0; f<nframes; f++) {
        //printf("read ADC\n");
        while(adc1_read(dummy_buf, nsam*2) == -1);  /* runs at Fs = 16kHz */

        //printf("Codec 2 enc\n");
	GPIOD->ODR = (1 << 13);
        codec2_encode(codec2, bits, pinbuf);
        pinbuf += nsam;
	GPIOD->ODR &= ~(1 << 13);
        //printf("Codec 2 dec\n");
	codec2_decode(codec2, outbuf, bits);

        //printf("write to DAC\n");
        while(dac1_write(dummy_buf, nsam*2) == -1); /* runs at Fs = 16kHz */
        //printf(".");
    }

    free(inbuf);
    free(outbuf);
    free(bits);
    codec2_destroy(codec2);
}

void gpio_init() {
    RCC->AHB1ENR |= RCC_AHB1ENR_GPIODEN; // enable the clock to GPIOD
    GPIOD->MODER = (1 << 26);            // set pin 13 to be general
                                         // purpose output
}

int main(int argc, char *argv[]) {
    SystemInit();
    gpio_init();
    machdep_profile_init ();
    adc_open(4*DAC_BUF_SZ);
    dac_open(4*DAC_BUF_SZ);

    printf("Starting power_ut\n");

    c2speedtest(CODEC2_MODE_1600, "stm_in.raw");

    printf("Finished\n");

    return 0;
}

