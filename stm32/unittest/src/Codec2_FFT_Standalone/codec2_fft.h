/*
 * codec2_fft.h
 *
 *  Created on: 17.09.2016
 *      Author: danilo
 */

#ifndef DRIVERS_FREEDV_CODEC2_FFT_H_
#define DRIVERS_FREEDV_CODEC2_FFT_H_

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#ifdef FDV_ARM_MATH
  #include "fdv_arm_math.h"
#else
    #define USE_KISS_FFT
#endif

#include "defines.h"
#include "comp.h"


typedef COMP    codec2_fft_cpx;

#ifdef USE_KISS_FFT
    #include "kiss_fft.h"
    typedef kiss_fft_cfg codec2_fft_cfg;
    typedef kiss_fft_scalar codec2_fft_scalar;
#else
  typedef float32_t codec2_fft_scalar;

  typedef struct {
        const arm_cfft_instance_f32* instance;
        int inverse;
    } codec2_fft_struct;
  typedef codec2_fft_struct* codec2_fft_cfg;
#endif



codec2_fft_cfg codec2_fft_alloc(int nfft, int inverse_fft, void* mem, size_t* lenmem);
void codec2_fft_free(codec2_fft_cfg cfg);


static inline void codec2_fft(codec2_fft_cfg cfg, codec2_fft_cpx* in, codec2_fft_cpx* out)
{

#ifdef USE_KISS_FFT
      kiss_fft(cfg, (kiss_fft_cpx*)in, (kiss_fft_cpx*)out);
#else
    memcpy(out,in,cfg->instance->fftLen*2*sizeof(float));
    arm_cfft_f32(cfg->instance,(float*)out,cfg->inverse,0);
    // TODO: this is not nice, but for now required to keep changes minimal
    // however, since main goal is to reduce the memory usage
    // we should convert to an in place interface
    // on PC like platforms the overhead of using the "inplace" kiss_fft calls
    // is neglectable compared to the gain in memory usage on STM32 platforms
    if (cfg->inverse)
    {
        arm_scale_f32((float*)out,cfg->instance->fftLen,(float*)out,cfg->instance->fftLen*2);
    }
#endif
}

void codec2_fft_inplace(codec2_fft_cfg cfg, codec2_fft_cpx* inout);


#endif
