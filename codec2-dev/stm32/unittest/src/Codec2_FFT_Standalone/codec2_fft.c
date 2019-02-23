/*
 * codec2_fft.c
 *
 *  Created on: 24.09.2016
 *      Author: danilo
 */


#define MALLOC malloc
#define CALLOC calloc
#define FREE free

#include "codec2_fft.h"

#ifdef USE_KISS_FFT
#include "_kiss_fft_guts.h"

#endif

void codec2_fft_free(codec2_fft_cfg cfg)
{
#ifdef USE_KISS_FFT
    KISS_FFT_FREE(cfg);
#else
    FREE(cfg);
#endif
}

codec2_fft_cfg codec2_fft_alloc(int nfft, int inverse_fft, void* mem, size_t* lenmem)
{
    codec2_fft_cfg retval;
#ifdef USE_KISS_FFT
    retval = kiss_fft_alloc(nfft, inverse_fft, mem, lenmem);
#else
    retval = MALLOC(sizeof(codec2_fft_struct));
    retval->inverse  = inverse_fft;
    switch(nfft)
    {
    case 128:
        retval->instance = &arm_cfft_sR_f32_len128;
        break;
    case 256:
        retval->instance = &arm_cfft_sR_f32_len256;
        break;
    case 512:
        retval->instance = &arm_cfft_sR_f32_len512;
        break;
//    case 1024:
//        retval->instance = &arm_cfft_sR_f32_len1024;
//        break;
    default:
        abort();
    }
    // retval->instance = arm_fft_cache_get(retval->instance);
#endif
    return retval;
}


// there is a little overhead for inplace kiss_fft but this is
// on the powerful platforms like the Raspberry or even x86 PC based ones
// not noticeable
// the reduced usage of RAM and increased performance on STM32 platforms
// should be worth it.
void codec2_fft_inplace(codec2_fft_cfg cfg, codec2_fft_cpx* inout)
{

#ifdef USE_KISS_FFT
    kiss_fft_cpx in[512];
    // decide whether to use the local stack based buffer for in
    // or to allow kiss_fft to allocate RAM
    // second part is just to play safe since first method
    // is much faster and uses less RAM
    if (cfg->nfft <= 512)
    {
        memcpy(in,inout,cfg->nfft*sizeof(kiss_fft_cpx));
        kiss_fft(cfg, in, (kiss_fft_cpx*)inout);
    }
    else
    {
        kiss_fft(cfg, (kiss_fft_cpx*)inout, (kiss_fft_cpx*)inout);
    }
#else
    arm_cfft_f32(cfg->instance,(float*)inout,cfg->inverse,1);
    if (cfg->inverse)
    {
        arm_scale_f32((float*)inout,cfg->instance->fftLen,(float*)inout,cfg->instance->fftLen*2);
    }

#endif
}
