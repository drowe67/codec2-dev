/*---------------------------------------------------------------------------*\

  FILE........: iir_tuner.c
  AUTHOR......: David Rowe
  DATE CREATED: 20 Feb 2015

  Filter/decimator function, broken out to this file so we can unit
  test easily.  

  Unit testing:
  
    ~/codec2-dev/stm32$ gcc -D__UNITTEST__ -Iinc src/iir_tuner.c -o iir_tuner -lm -Wal
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

/* Filter coefficients of IIR tuner (BETA1) and FIR equaliser (BETA2).
   Note neat trick to relate BETA2 to BETA1 by the decimation rate */

#define BETA1                    0.999
#define BETA2                    (1.0 - (1.0-BETA1)*ADC_TUNER_M)

/* filter states - we keep them global due to the need for speed */

float y_2, y_1, z_2, z_1;

/*
   ADC -> signed conversion - IIR BPF - Decimate - FIR Equaliser -> FIFO
*/

void inline iir_tuner(float dec_buf[], unsigned short adc_buf[]) {
    int i, j, k;
    float x, y, z;

    for(i=0, j=0; i<ADC_TUNER_BUF_SZ/2; j++) {

        /* IIR BPF centred at Fs/4.  All your MIPs are belong to this
           loop. */

        for(k=0; k<ADC_TUNER_M; k++,i++) {
            x = (int)adc_buf[i] - 32768;
            y = x - BETA1*y_2;
            y_2 = y_1;
            y_1 = y;
        }

        /* Equaliser FIR filter, notch at Fs/(4*ADC_TUNER_M) to smooth out 
           IIR BPF passband response */

        z = y + BETA2*z_2;
        dec_buf[j] = z;
        z_2 = z_1;
        z_1 = y;
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
    float          dec_s[NOUT];
    FILE          *f;
    int            i,j;

    f1 = 500E3;
    f2 = f1 + 8E3;       /* wanted */
    f3 = f1 - 7E3;       /* wanted */
    f4 = f1 - 207E3;     /* out of band, should be greatly attenuated */

    for(i=0; i<NIN; i++)
        s[i] = 32767;
    synth_line(s, f2, 0.1, NIN);
    synth_line(s, f3, 0.05, NIN);
    synth_line(s, f4, 0.1, NIN);
    for(i=0, j=0; i<NIN; i+=ADC_TUNER_BUF_SZ/2, j+=ADC_TUNER_N/2) {
        iir_tuner(&dec_s[j], &s[i]);
    }
    
    f = fopen("iir_tuner_s.txt", "wt");  assert(f != NULL);
    for(i=0; i<NIN; i++)
        fprintf(f, "%d\n", s[i]);
    fprintf(f, "\n");
    fclose(f);

    f = fopen("iir_tuner.txt", "wt");  assert(f != NULL);
    for(i=0; i<NOUT; i++)
        fprintf(f, "%f\n", dec_s[i]);
    fprintf(f, "\n");
    fclose(f);

    return 0;
}

#endif
