#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>

#include "comp.h"
#include "codec2_fft.h"
#include "kiss_fft.h"

#define FFT_ENC  512
#define M_PITCH_S  0.0400 

int main(int argc, char *argv[]) {

    int           m_pitch;
    int           nw;
    float        *w;
    COMP          W[FFT_ENC];
    codec2_fft_cfg  fft_fwd_cfg;
    int i, j;

    m_pitch = floor(8000*M_PITCH_S);
    nw = 279;

    w = (float*)malloc(m_pitch*sizeof(float));

#ifdef USE_KISS_FFT
    fft_fwd_cfg  = kiss_fft_alloc(FFT_ENC, 0, NULL, NULL);
#else
    fft_fwd_cfg  = malloc(sizeof(codec2_fft_struct));
    fft_fwd_cfg ->inverse  = 0;
    fft_fwd_cfg ->instance = &arm_cfft_sR_f32_len512;
#endif

    COMP  wshift[FFT_ENC];

    float m;
    m = 0.0;
    for(i=0; i<m_pitch/2-nw/2; i++)
      w[i] = 0.0;
    for(i=m_pitch/2-nw/2,j=0; i<m_pitch/2+nw/2; i++,j++) {
      w[i] = 0.5 - 0.5*cosf(TWO_PI*j/(nw-1));
      m += w[i]*w[i];
    }
    for(i=m_pitch/2+nw/2; i<m_pitch; i++)
      w[i] = 0.0;

    m = 1.0/sqrtf(m*FFT_ENC);
    for(i=0; i<m_pitch; i++) {
      w[i] *= m;
    }


    for(i=0; i<FFT_ENC; i++) {
      wshift[i].real = 0.0;
      wshift[i].imag = 0.0;
    }
    for(i=0; i<nw/2; i++)
      wshift[i].real = w[i+m_pitch/2];
    for(i=FFT_ENC-nw/2,j=m_pitch/2-nw/2; i<FFT_ENC; i++,j++)
     wshift[i].real = w[j];


    printf("static const COMP wshift[] = {\n    ");
    j = 0;
    for (i=0; i<FFT_ENC; i++) {
        if (j>0) printf(" ");
        printf("{%f,%f}", (double)wshift[i].real, (double)wshift[i].imag);
        if (i<(FFT_ENC-1)) printf(",");
        if (++j >=4) {
            printf("\n    ");
            j = 0;
            }
        }
    printf("};\n\n");

    /////////////////
    kiss_fft(fft_fwd_cfg, (kiss_fft_cpx*)wshift, (kiss_fft_cpx*)W);
    /////////////////


    printf("static const float W_expect[] = {\n    ");
    j = 0;
    for (i=0; i<FFT_ENC; i++) {
        if (j>0) printf(" ");
        printf("%f", (double)W[i].real);
        if (i<(FFT_ENC-1)) printf(",");
        if (++j >=8) {
            printf("\n    ");
            j = 0;
            }
        }
    printf("};\n\n");

    fclose(stdout);
    fclose(stderr);

    return(0);
}

/* vi:set ts=4 et sts=4: */
