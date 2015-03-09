 /*---------------------------------------------------------------------------*\

  FILE........: iir_duc.c
  AUTHOR......: Brady O'Brien
  DATE CREATED: 6 Mar 2015

  Interapolator/Filter for IF upconversion

  Unit testing:
  
    ~/codec2-dev/stm32$ gcc -D__UNITTEST__ -Iinc src/iir_duc.c -o iir_duc -lm -Wall
    ~/codec2-dev/stm32$ ./iir_duc

\*---------------------------------------------------------------------------*/

/*
  Copyright (C) 2015 Brady O'Brien

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

#include "stm32f4_dacduc.h"
#include "iir_duc.h"

#define BETA1                    0.99002			// B1MUL/(2**B1SHFT)
#define B1MUL			 32441	
#define B1SMUL			 -38328
#define B1SHFT			 15				// 10 bits gives us plenty of headroom between 31 bits of int and 14 bits of ADC
#define B2MUL			 24593				// This actually matches BETA2 exactly with the supplied BETA1
#define B2SHFT			 15				// 10 is also the lowest we can go without beta1=1
#define BETA2                    (1.0 - (1.0-BETA1)*DUC_M)	// B2MUL/(2**B2SHFT)
#define IN_SCALE                 2.0                            //Input scaling factor. Should be as large as the amplitude of the incoming samples
#define DAC_SCALE                4096                           //Maximum output to DAC
#define DAC_SCALE_2		 2040


//IIR and FIR filter states. Global for go fast.
float f_1,f_2,f;
int   n_1,n_2,n;

/*
   Upconvert and bandpass filter a chunk of spectrum from Fs/M to Fs. We're going for 700khz here.
   modin needs to be DUC_N long and dac_out needs to be DUC_N*DUC_M long. This 
*/

void iir_upconv(float modin[], unsigned short dac_out[]){
    int i,j,k;
    int m;
    k=0;
    //Iterate through input samples and apply pre-eq FIR, interpolate, and apply BPF IIR
    for(i=0;i<DUC_N;i++){
        f = modin[i]+f_2*BETA2;
        f_2 = f_1;
        f_1 = modin[i];                                             //Scale fir output and convert to fixed.
        m = (int)((f/(IN_SCALE))*DAC_SCALE_2);                      //Scale fir output and convert to fixed
        n = (m) + ((B1SMUL*n_1)>>B1SHFT) - ((B1MUL*n_2)>>B1SHFT);   //Apply one cycle of IIR. This feeds the fir-ed sample into the output filter
        n_2 = n_1;
        n_1 = n;
        dac_out[k]=(unsigned short)(n+DAC_SCALE_2);
        k++;
        //now do the rest of the filtering. Because we're zero-stuffing we can neglect the sample from the fir filter.
        for(j=1;j<DUC_M;j++,k++){
            n = ((B1SMUL*n_1)>>B1SHFT) - ((B1MUL*n_2)>>B1SHFT);
            n_2 = n_1;
            n_1 = n;
            dac_out[k]=(unsigned short)((n)+DAC_SCALE_2);
        }
    }
}


#ifdef __UNITTEST__

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

#define FS     80000
#define AMP_MAX 1

#define NOUT_BUFS    401
#define NIN          (NOUT_BUFS*DUC_N)
#define NOUT         (NIN*DUC_M)

void synth_line(float us[], float f, float amp, int n) {
    float w, sam;
    int   i;

    w = 2*M_PI*f/(float)FS;

    for(i=0; i<n; i++) {
        sam = amp*AMP_MAX*cos(w*i);
        us[i] += sam;
    }
}

float 	   s[NIN];
float      fout[NIN];
unsigned short todac[NOUT];

int main(void) {
    float          f1,f2,f3;
    FILE          *f;
    int            i;

    f1 = 20E3;	          /* center of passband */
    f2 = f1;	  /* wanted */
    f3 = f1 + 7E3;        /* wanted */


    for(i=0; i<NIN; i++)
        s[i] = 0.01;
    synth_line(s, f2, 0.5, NIN);
    synth_line(s, f3, 0.5, NIN);
    for(i=0;i<NOUT_BUFS;i++)
        iir_upconv(&s[i*(DUC_N)],&todac[i*(DUC_N*DUC_M)]);
    
    f = fopen("iir_duc_s.txt", "wt");  assert(f != NULL);
    for(i=DUC_N; i<NIN; i++)
        fprintf(f, "%f\n", s[i]);
    fprintf(f, "\n");
    fclose(f);

    f = fopen("iir_duc_f.txt", "wt");  assert(f != NULL);
    for(i=DUC_N; i<NIN; i++)
        fprintf(f, "%f\n", fout[i]);
    fprintf(f, "\n");
    fclose(f);

    f = fopen("iir_duc.txt", "wt");  assert(f != NULL);
    for(i=DUC_N*DUC_M; i<NOUT; i++)
        fprintf(f, "%d\n", todac[i]);
    fprintf(f, "\n");
    fclose(f);

    return 0;
}

#endif



