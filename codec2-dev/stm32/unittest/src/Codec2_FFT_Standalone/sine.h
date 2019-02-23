#ifndef __SINE__
#define __SINE__

#include "defines.h"
#include "comp.h"
#include "codec2_fft.h"

C2CONST c2const_create(int Fs, float framelength_ms);

void make_analysis_window(C2CONST *c2const, codec2_fft_cfg fft_fwd_cfg, float w[], COMP W[]);

#endif
