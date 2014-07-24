/* 
   t16_8.c
   David Rowe
   May 10 2012

   Unit test for 16 to 8 kHz sample rate conversion functions.  I
   evaluated output by plotting using Octave and looking for jaggies:

     pl("../unittest/out16.raw",1,3000)
     pl("../unittest/out8.raw",1,3000)

   Listening to it also shows up anything nasty:

     $ play -s -2 -r 16000 out16.raw
     $ play -s -2 -r 8000 out8.raw

  */

#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "codec2_fdmdv.h"

#define N8                        160 /* procssing buffer size at 8 kHz */
#define N16             (N8*FDMDV_OS)
#define MEM8 (FDMDV_OS_TAPS/FDMDV_OS)
#define FRAMES                     50
#define TWO_PI            6.283185307
#define FS                       8000

#define SINE

int main() {
    float in8k[MEM8 + N8];
    float out16k[N16];
    short out16k_short[N16];
    FILE *f16;

    float in16k[FDMDV_OS_TAPS + N16];
    float out8k[N16];
    short out8k_short[N8];
    FILE *f8;

    int i,f,t,t1;
    float freq = 800.0;

    f16 = fopen("out16.raw", "wb");
    assert(f16 != NULL);
    f8 = fopen("out8.raw", "wb");
    assert(f8 != NULL);
    
    /* clear filter memories */

    for(i=0; i<MEM8; i++)
	in8k[i] = 0.0;
    for(i=0; i<FDMDV_OS_TAPS; i++)
	in16k[i] = 0.0;

    t = t1 = 0;
    for(f=0; f<FRAMES; f++) {

#ifdef DC
	for(i=0; i<N8; i++)
	    in8k[MEM8+i] = 16000.0;
#endif
#ifdef SINE
	for(i=0; i<N8; i++,t++)
	    in8k[MEM8+i] = 16000.0*cos(TWO_PI*t*freq/FS);
#endif

	/* upsample  */

	fdmdv_8_to_16(out16k, &in8k[MEM8], N8);
	/*
	for(i=0; i<MEM8; i++)
	    in8k[i] = in8k[i+N8];
	*/

	/* save 16k to disk for plotting and check out */

	for(i=0; i<N16; i++)
	    out16k_short[i] = (short)out16k[i];
	fwrite(out16k_short, sizeof(short), N16, f16);
	
	/* add a 10 kHz spurious signal for fun, we want down sampler to
	   knock this out */

	for(i=0; i<N16; i++,t1++)
	    in16k[i+FDMDV_OS_TAPS] = out16k[i] + 16000.0*cos(TWO_PI*t1*1E4/FS);

	/* downsample */

	fdmdv_16_to_8(out8k, &in16k[FDMDV_OS_TAPS], N8);
	/*
	for(i=0; i<FDMDV_OS_TAPS; i++)
	    in16k[i] = in16k[i+N16];
	*/

	/* save 8k to disk for plotting and check out */

	for(i=0; i<N8; i++)
	    out8k_short[i] = (short)out8k[i];
	fwrite(out8k_short, sizeof(short), N8, f8);
	
    }

    fclose(f16);
    fclose(f8);
    return 0;

}
