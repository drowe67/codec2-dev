/*
   t24_8.c
   Mooneer Salem
   16 Feb 2023

   Unit test for 24 to 8 kHz sample rate conversion functions.
  */

#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "codec2_fdmdv.h"

#define N8         180        /* processing buffer size at 8 kHz */
#define N24        (N8*FDMDV_OS_24)
#define MEM8       FDMDV_OS_TAPS_24_8K
#define FRAMES     50
#define TWO_PI     6.283185307
#define FS         48000

int main() {
    short in8k[MEM8 + N8];
    float out24k[N24];
    FILE *f24;

    float in24k[FDMDV_OS_TAPS_24K + N24];
    short out8k[N8];
    FILE *f8, *f8in;

    int i,f,t,t1;
    float freq = 800.0;

    f24 = fopen("out24.raw", "wb");
    assert(f24 != NULL);
    f8 = fopen("out8.raw", "wb");
    assert(f8 != NULL);
    f8in = fopen("in8.raw", "wb");
    assert(f8in != NULL);

    /* clear filter memories */

    for(i=0; i<MEM8; i++)
        in8k[i] = 0.0;
    for(i=0; i<FDMDV_OS_TAPS_24K; i++)
        in24k[i] = 0.0;

    t = t1 = 0;
    for(f=0; f<FRAMES; f++) {

        for(i=0; i<N8; i++,t++)
	        in8k[MEM8+i] = 16000.0*cos(TWO_PI*t*freq/(FS/FDMDV_OS_24));

        fwrite(in8k, sizeof(short), N8, f8in);

        /* upsample  */

        fdmdv_8_to_24(out24k, &in8k[MEM8], N8);

        /* save 24k to disk for plotting and check out */
        fwrite(out24k, sizeof(float), N24, f24);

    	/* add a 10 kHz spurious signal for fun, we want down sampler to
    	   knock this out */
        for(i=0; i<N24; i++,t1++)
	        in24k[i+FDMDV_OS_TAPS_24K] = out24k[i] + FDMDV_SHORT_TO_FLOAT * 16000.0*cos(TWO_PI*t1*1E4/FS);

	    /* downsample */
        fdmdv_24_to_8(out8k, &in24k[FDMDV_OS_TAPS_24K], N8);

        /* save 8k to disk for plotting and check out */
        fwrite(out8k, sizeof(short), N8, f8);

    }

    fclose(f24);
    fclose(f8);
    fclose(f8in);
    return 0;

}
