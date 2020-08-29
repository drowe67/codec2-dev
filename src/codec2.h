/*---------------------------------------------------------------------------*\

  FILE........: codec2.h
  AUTHOR......: David Rowe
  DATE CREATED: 21 August 2010

  Codec 2 fully quantised encoder and decoder functions.  If you want use
  Codec 2, these are the functions you need to call.

\*---------------------------------------------------------------------------*/

/*
  Copyright (C) 2010 David Rowe

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

#ifndef __CODEC2__
#define  __CODEC2__

#include <codec2/version.h>

#ifdef __cplusplus
  extern "C" {
#endif

#if defined _WIN32 || defined __CYGWIN__
#ifdef __GNUC__
#define API_EXPORT __attribute__ ((dllexport))
#else
#define API_EXPORT __declspec(dllexport)
#endif
#else
#define API_EXPORT __attribute__ ((visibility ("default")))
#endif

#define CODEC2_MODE_3200 	0
#define CODEC2_MODE_2400 	1
#define CODEC2_MODE_1600 	2
#define CODEC2_MODE_1400 	3
#define CODEC2_MODE_1300 	4
#define CODEC2_MODE_1200 	5
#define CODEC2_MODE_700C 	8
#define CODEC2_MODE_450 	10
#define CODEC2_MODE_450PWB 	11

#ifndef CODEC2_MODE_EN_DEFAULT
#define CODEC2_MODE_EN_DEFAULT 1
#endif

// by default we enable all modes
// disable during compile time with -DCODEC2_MODE_1600_EN=0
// all but CODEC2 1600 are enabled then

//or the other way round
// -DCODEC2_MODE_EN_DEFAULT=0 -DCODEC2_MODE_1600_EN=1
// only CODEC2 Mode 1600

#if !defined(CODEC2_MODE_3200_EN)
        #define CODEC2_MODE_3200_EN CODEC2_MODE_EN_DEFAULT
#endif
#if !defined(CODEC2_MODE_2400_EN)
        #define CODEC2_MODE_2400_EN CODEC2_MODE_EN_DEFAULT
#endif
#if !defined(CODEC2_MODE_1600_EN)
        #define CODEC2_MODE_1600_EN CODEC2_MODE_EN_DEFAULT
#endif
#if !defined(CODEC2_MODE_1400_EN)
        #define CODEC2_MODE_1400_EN CODEC2_MODE_EN_DEFAULT
#endif
#if !defined(CODEC2_MODE_1300_EN)
        #define CODEC2_MODE_1300_EN CODEC2_MODE_EN_DEFAULT
#endif
#if !defined(CODEC2_MODE_1200_EN)
        #define CODEC2_MODE_1200_EN CODEC2_MODE_EN_DEFAULT
#endif
#if !defined(CODEC2_MODE_700C_EN)
        #define CODEC2_MODE_700C_EN CODEC2_MODE_EN_DEFAULT
#endif
#if !defined(CODEC2_MODE_450_EN)
        #define CODEC2_MODE_450_EN CODEC2_MODE_EN_DEFAULT
#endif
#if !defined(CODEC2_MODE_450PWB_EN)
        #define CODEC2_MODE_450PWB_EN CODEC2_MODE_EN_DEFAULT
#endif

#define CODEC2_MODE_ACTIVE(mode_name, var)  ((mode_name##_EN) == 0 ? 0: (var) == mode_name)

struct CODEC2;

API_EXPORT struct CODEC2 *  codec2_create(int mode);
API_EXPORT void codec2_destroy(struct CODEC2 *codec2_state);
API_EXPORT void codec2_encode(struct CODEC2 *codec2_state, unsigned char * bits, short speech_in[]);
API_EXPORT void codec2_decode(struct CODEC2 *codec2_state, short speech_out[], const unsigned char *bits);
API_EXPORT void codec2_decode_ber(struct CODEC2 *codec2_state, short speech_out[], const unsigned char *bits, float ber_est);
API_EXPORT int  codec2_samples_per_frame(struct CODEC2 *codec2_state);
API_EXPORT int  codec2_bits_per_frame(struct CODEC2 *codec2_state);

API_EXPORT void codec2_set_lpc_post_filter(struct CODEC2 *codec2_state, int enable, int bass_boost, float beta, float gamma);
API_EXPORT int  codec2_get_spare_bit_index(struct CODEC2 *codec2_state);
API_EXPORT int  codec2_rebuild_spare_bit(struct CODEC2 *codec2_state, char unpacked_bits[]);
API_EXPORT void codec2_set_natural_or_gray(struct CODEC2 *codec2_state, int gray);
API_EXPORT void codec2_set_softdec(struct CODEC2 *c2, float *softdec);
API_EXPORT float codec2_get_energy(struct CODEC2 *codec2_state, const unsigned char *bits);
      
// support for ML and VQ experiments
API_EXPORT void codec2_open_mlfeat(struct CODEC2 *codec2_state, char *feat_filename, char *model_filename);
API_EXPORT void codec2_load_codebook(struct CODEC2 *codec2_state, int num, char *filename);
API_EXPORT float codec2_get_var(struct CODEC2 *codec2_state);
API_EXPORT float *codec2_enable_user_ratek(struct CODEC2 *codec2_state, int *K);

// 700C post filter and equaliser
API_EXPORT void codec2_700c_post_filter(struct CODEC2 *codec2_state, int en);
API_EXPORT void codec2_700c_eq(struct CODEC2 *codec2_state, int en);
      
#ifdef __cplusplus
}
#endif

#endif

