/* 
   t16__short_8.c
   David Rowe
   19 August 2014

   Unit test for 16 to 8 kHz sample rate conversion functions.  I
   evaluated output by plotting using Octave and looking for jaggies:

     pl("../unittest/out16_short.raw",1,3000)
     pl("../unittest/out8.raw",1,3000)

   Listening to it also shows up anything nasty:

     $ play -s -2 -r 16000 out16_short.raw
     $ play -s -2 -r 8000 out8.raw

  */

#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "codec2_fdmdv.h"

#define N8                        160 /* procssing buffer size at 8 kHz */
#define N16             (N8*FDMDV_OS)
#define FRAMES                     100
#define TWO_PI            6.283185307
#define FS                      16000

#define SINE

int main() {
    float in8k[FDMDV_OS_TAPS_8K + N8];
    short out16k_short[N16];
    FILE *f16;

    short in16k_short[FDMDV_OS_TAPS_16K + N16];
    float out8k[N16];
    short out8k_short[N8];
    FILE *f8;

    int i,f,t,t1;
    float freq = 800.0;

    f16 = fopen("out16_short.raw", "wb");
    assert(f16 != NULL);
    f8 = fopen("out8.raw", "wb");
    assert(f8 != NULL);
    
    /* clear filter memories */

    for(i=0; i<FDMDV_OS_TAPS_8K; i++)
	in8k[i] = 0.0;
    for(i=0; i<FDMDV_OS_TAPS_16K; i++)
	in16k_short[i] = 0;

    t = t1 = 0;
    for(f=0; f<FRAMES; f++) {

#ifdef DC
	for(i=0; i<N8; i++)
	    in8k[FDMDV_OS_TAPS_8K+i] = 16000.0;
#endif
#ifdef SINE
	for(i=0; i<N8; i++,t++)
	    in8k[FDMDV_OS_TAPS_8K+i] = 8000.0*cos(TWO_PI*t*freq/FS);
#endif

	/* upsample  */

	fdmdv_8_to_16_short(out16k_short, &in8k[FDMDV_OS_TAPS_8K], N8);

	fwrite(out16k_short, sizeof(short), N16, f16);
	
	/* add a 6 kHz spurious signal for fun, we want down sampler to
	   knock this out */

	for(i=0; i<N16; i++,t1++)
	    in16k_short[i+FDMDV_OS_TAPS_16K] = out16k_short[i] + 8000.0*cos(TWO_PI*t1*6000.0/FS);

	/* downsample */

	fdmdv_16_short_to_8(out8k, &in16k_short[FDMDV_OS_TAPS_16K], N8);

	/* save 8k to disk for plotting and check out */

	for(i=0; i<N8; i++)
	    out8k_short[i] = (short)out8k[i];
	fwrite(out8k_short, sizeof(short), N8, f8);
	
    }

    fclose(f16);
    fclose(f8);
    return 0;

}
