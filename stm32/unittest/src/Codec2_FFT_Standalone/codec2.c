
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "defines.h"
#include "codec2.h"
#include "codec2_internal.h"
#include "codec2_fft.h"
#include "sine.h"


#define MALLOC malloc
#define CALLOC calloc
#define FREE free


struct CODEC2 * codec2_create(int mode)
{
    struct CODEC2 *c2;

    c2 = (struct CODEC2*)MALLOC(sizeof(struct CODEC2));
    if (c2 == NULL)
	return NULL;

    c2->mode = mode;

    /* store constants in a few places for convenience */
    
    c2->c2const = c2const_create(8000, N_S);
    c2->Fs = c2->c2const.Fs;
    int m_pitch = c2->m_pitch = c2->c2const.m_pitch;

    c2->w = (float*)MALLOC(m_pitch*sizeof(float));
    if (c2->w == NULL) {
        FREE(c2->Pn);
        FREE(c2->Sn_);
	return NULL;
    }

    c2->fft_fwd_cfg = codec2_fft_alloc(FFT_ENC, 0, NULL, NULL);

    make_analysis_window(&c2->c2const, c2->fft_fwd_cfg, c2->w,c2->W);

    return c2;
}
