 /*---------------------------------------------------------------------------*\

  FILE........: iir_tuner.c
  AUTHOR......: David Rowe
  DATE CREATED: 20 Feb 2015

  Filter/decimator function, broken out to this file so we can unit
  test easily.

  Unit testing:

    ~/codec2-dev/stm32$ gcc -D__UNITTEST__ -Iinc src/iir_tuner.c -o iir_tuner -lm -Wall
    ~/codec2-dev/stm32$ ./iir_tuner

\*---------------------------------------------------------------------------*/

/*
  Copyright (C) 2015 David Rowe

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

#ifdef __UNITTEST__

#include <assert.h>
#include <math.h>
#include <stdio.h>

#endif

#include "stm32f4_adc_tuner.h"
#include "iir_tuner.h"

/* Filter coefficients of IIR tuner (BETA1) and FIR equaliser (BETA2).
   Note neat trick to relate BETA2 to BETA1 by the decimation rate */

#define BETA1                    .9990234375			// B1MUL/(2**B1SHFT)
#define B1MUL			 1023
#define B1SMUL			 1204
#define B1SHFT			 10				// 10 bits gives us plenty of headroom between 31 bits of int and 14 bits of ADC
#define B2MUL			 979				// This actually matches BETA2 exactly with the supplied BETA1
#define B2SHFT			 10				// 10 is also the lowest we can go without beta1=1
#define BETA2                    (1.0 - (1.0-BETA1)*ADC_TUNER_M)// B2MUL/(2**B2SHFT)

#define FIXED_IIR                                               //Define this to compile a fixed point IIR filter

/* filter states - we keep them global due to the need for speed */

#ifdef FIXED_IIR
int n_2, n_1, o_2, o_1;
#else
float y_2, y_1, z_2, z_1;
#endif

/*
   ADC -> signed conversion - IIR BPF - Decimate - FIR Equaliser -> FIFO
*/

void iir_tuner(
               float          dec_50[],   // ADC_TUNER_N/2 output samples
               unsigned short adc_buf[]   // ADC_TUNER_BUF_SZ/2 input samples
) 
{
    int i, j, k;
#ifndef FIXED_IIR
    float x, y, z;
#endif
    int n, m, o;

    for(i=0, j=0; i<ADC_TUNER_BUF_SZ/2; j++) {

        /* IIR BPF centred at Fs/4.  All your MIPs are belong to this
           loop. */
        for(k=0; k<ADC_TUNER_M; k++,i++) {
            #ifdef FIXED_IIR
            m = (int)adc_buf[i];
            n = m - ((B1SMUL*n_1)>>B1SHFT) - ((B1MUL*n_2)>>B1SHFT);
            n_2 = n_1;
            n_1 = n;
            #else
	    x = (int)adc_buf[i] - 32768;
	    y = x - (BETA1*y_2);
	    y_2 = y_1;
            y_1 = y;
            #endif
        }

        /* Equaliser FIR filter, notch at Fs/(4*ADC_TUNER_M) to smooth out
           IIR BPF passband response */
        #ifdef FIXED_IIR
	o = n + ((B2MUL*o_2)>>B2SHFT);
	dec_50[j] = (float)o;
	o_2 = o_1;
	o_1 = n;
        #else
	z=y+BETA2*z_2;
	dec_50[j] = z;
        z_2 = z_1;
        z_1 = y;
        #endif

    }
}


/* BPF at 12.5 kHz +/- 2000 Hz, and decimate down to Fs = 10kHz */

static float fir_50_to_10[];
void iir_tuner_dec_50_to_10(float dec_10[], float dec_50[], int n) {
    int   i,j,k;
    float acc;

    for(i=0,k=0; i<n; i+=5,k++) {
        acc = 0.0;
        for(j=0; j<IIR_TUNER_DEC_50_10_FILT_MEM; j++)
            acc += dec_50[i-j]*fir_50_to_10[j];
        dec_10[k] = acc;
    }

}


#ifdef __UNITTEST__

#define FS      2000000
#define AMP_MAX 32767

#define NOUT_BUFS    100
#define NOUT         (NOUT_BUFS*ADC_TUNER_N)
#define NIN          (NOUT*ADC_TUNER_M)

void synth_line(unsigned short us[], float f, float amp, int n) {
    float w, sam;
    int   i;

    w = 2*M_PI*f/(float)FS;

    for(i=0; i<n; i++) {
        sam = amp*AMP_MAX*cos(w*i);
        us[i] += (unsigned short)(sam + 0.5);
    }
}


int main(void) {
    float          f1,f2,f3,f4;
    unsigned short s[NIN];
    float          dec_50[IIR_TUNER_DEC_50_10_FILT_MEM+NOUT];
    float          dec_10[NOUT/5];
    FILE          *f;
    int            i,j,k;
    short          dec_10_short;

    /* test Fs=2E6 unsigned short to Fs=50E3 float tuner/resampler -----------------------*/

    f1 = 700E3;
    f2 = f1 + 8E3;       /* wanted */
    f3 = f1 - 7E3;       /* wanted */
    f4 = f1 - 207E3;     /* out of band, should be greatly attenuated */

    for(i=0; i<NIN; i++)
        s[i] = 32767;
    synth_line(s, f1, 1, NIN);
    //synth_line(s, f3, 0.1, NIN);
    //synth_line(s, f4, 0.2, NIN);
    for(i=0, j=0; i<NIN; i+=ADC_TUNER_BUF_SZ/2, j+=ADC_TUNER_N/2) {
        iir_tuner(&dec_50[j], &s[i]);
    }

    f = fopen("iir_tuner_s.txt", "wt");  assert(f != NULL);
    for(i=0; i<NIN; i++)
        fprintf(f, "%d\n", s[i]);
    fprintf(f, "\n");
    fclose(f);

    f = fopen("iir_tuner.txt", "wt");  assert(f != NULL);
    for(i=0; i<NOUT; i++)
        fprintf(f, "%f\n", dec_50[i]);
    fprintf(f, "\n");
    fclose(f);

    /* test FS=2E6 unsigned short -> Fs=10kHz short ---------------------------------------------*/

    for(i=0; i<NIN; i++)
        s[i] = 32767;
    for(i=1; i<NIN; i+=4)
        s[i] += 32767;
    for(i=3; i<NIN; i+=4)
        s[i] -= 32767;

    for(i=0, j=0, k=0; i<NIN; i+=ADC_TUNER_BUF_SZ/2, j+=ADC_TUNER_N/2, k+=(ADC_TUNER_N/2)/5) {
        iir_tuner(&dec_50[IIR_TUNER_DEC_50_10_FILT_MEM+j], &s[i]);
        iir_tuner_dec_50_to_10(&dec_10[k], &dec_50[IIR_TUNER_DEC_50_10_FILT_MEM+j], ADC_TUNER_N/2);
    }

    f = fopen("iir_tuner2.txt", "wt");  assert(f != NULL);
    for(i=0; i<NOUT/5; i++) {
        dec_10_short = dec_10[i]/ADC_TUNER_M;
        fprintf(f, "%d\n", dec_10_short);
    }
    fprintf(f, "\n");
    fclose(f);

    return 0;
}

#endif

/* Band pass FIR filter coefficents centred on Fs/4, used before Fs=50kHz to Fs=10kHz */

static float fir_50_to_10[] = {
    -1.71502876e-07,
    -3.93029078e-05,
    -5.30743362e-04,
    1.17938704e-04,
    1.09727519e-03,
    -1.90605585e-04,
    -1.61350037e-03,
    2.37746793e-04,
    1.86947117e-03,
    -2.28459776e-04,
    -1.56457257e-03,
    1.33627883e-04,
    4.28669971e-04,
    5.51555269e-05,
    1.60217953e-03,
    -3.09989338e-04,
    -4.22595074e-03,
    5.63604865e-04,
    6.71504730e-03,
    -7.23797964e-04,
    -8.02135199e-03,
    7.03165104e-04,
    7.05924883e-03,
    -4.54750532e-04,
    -3.11212393e-03,
    2.24463518e-08,
    -3.75334414e-03,
    5.63496992e-04,
    1.24113249e-02,
    -1.08145162e-03,
    -2.06694342e-02,
    1.38572694e-03,
    2.55955103e-02,
    -1.34897285e-03,
    -2.41472078e-02,
    9.32473244e-04,
    1.39715800e-02,
    -2.09763900e-04,
    5.83111135e-03,
    -6.44614872e-04,
    -3.42028021e-02,
    1.40049434e-03,
    6.80569757e-02,
    -1.83898122e-03,
    -1.02710848e-01,
    1.82014182e-03,
    1.32754277e-01,
    -1.32936318e-03,
    -1.53163914e-01,
    4.86473969e-04,
    1.60393264e-01,
    4.86473969e-04,
    -1.53163914e-01,
    -1.32936318e-03,
    1.32754277e-01,
    1.82014182e-03,
    -1.02710848e-01,
    -1.83898122e-03,
    6.80569757e-02,
    1.40049434e-03,
    -3.42028021e-02,
    -6.44614872e-04,
    5.83111135e-03,
    -2.09763900e-04,
    1.39715800e-02,
    9.32473244e-04,
    -2.41472078e-02,
    -1.34897285e-03,
    2.55955103e-02,
    1.38572694e-03,
    -2.06694342e-02,
    -1.08145162e-03,
    1.24113249e-02,
    5.63496992e-04,
    -3.75334414e-03,
    2.24463518e-08,
    -3.11212393e-03,
    -4.54750532e-04,
    7.05924883e-03,
    7.03165104e-04,
    -8.02135199e-03,
    -7.23797964e-04,
    6.71504730e-03,
    5.63604865e-04,
    -4.22595074e-03,
    -3.09989338e-04,
    1.60217953e-03,
    5.51555269e-05,
    4.28669971e-04,
    1.33627883e-04,
    -1.56457257e-03,
    -2.28459776e-04,
    1.86947117e-03,
    2.37746793e-04,
    -1.61350037e-03,
    -1.90605585e-04,
    1.09727519e-03,
    1.17938704e-04,
    -5.30743362e-04,
    -3.93029078e-05,
    -1.71502876e-07
};

