/*---------------------------------------------------------------------------*\

  FILE........: ofdm_internal.h
  AUTHORS.....: David Rowe & Steve Sampson
  DATE CREATED: June 2017

  OFDM Internal definitions.

\*---------------------------------------------------------------------------*/

/*
  Copyright (C) 2017 David Rowe

  All rights reserved.

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU Lesser General Public License version 2, as
  published by the Free Software Foundation.  This program is
  distributed in the hope that it will be useful, but WITHOUT ANY
  WARRANTY; without even the implied warranty of MERCHANTABILITY or
  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public
  License for more details.

  You should have received a copy of the GNU Lesser General Public License
  along with this program; if not, see <http://www.gnu.org/licenses/>.
*/

#ifndef OFDM_INTERNAL_H
#define OFDM_INTERNAL_H

#ifdef __cplusplus
extern "C" {
#endif

#include <complex.h>
#include <stdbool.h>
#include <stdint.h>

#include "codec2_ofdm.h"
#include "kiss_fft.h"

#ifndef M_PI
#define M_PI        3.14159265358979323846f  /* math constant */
#endif

#define TAU         (2.0f * M_PI)            /* mathematical constant */

struct OFDM_CONFIG{
    int32_t            Nc;                /* N carriers*/
    float              Ts;                /* Symbol time (S) (no CP) */
    float              Rs;                /* Symbol rate (S) (no CP) */
    int32_t            Fs;                /* Sample rate */
    int32_t            bps;               /* Bits per symbol */
    int32_t            Tcp;               /* Cyclic Prefix time (S) */
    int32_t            Ns;                /* Symbols per frame, including pilot */
    int32_t            Fcenter;           /* Center frequency */
    int32_t            M;                 /* Samples per symbol (no CP) */
    int32_t            Ncp;               /* Samples of cyclic prefix */
    int32_t            FtWindowWidth;     /* ? */
    int32_t            BitsPerFrame;      /* Bits in a frame */
    int32_t            SampsPerFrame;     /* Samples in a frame */
    int32_t            SampsPerFrameMax;  /* Maximum samples per demod cycle */
    int32_t            RxBufSize;         /* RX buffer size */
    int32_t            RowsPerFrame;      /* Symbol depth per frame, no prefix */
};

struct OFDM {
    struct OFDM_CONFIG config;            /* OFDM configuration (see above) */
    float foff_est_gain;                  /* ? */
    float foff_est_hz;                    /* Estimated frequency offset */

    int verbose;                          /* Printing verbosity level */
    int sample_point;                     /* Fine timing offset */
    int timing_est;                       /* Ditto (?) (Coarse timing, maybe?) */
    int nin;                              /* Samples expected next demod cycle */
    int frame_point;
    int sync_count;

    bool timing_en;                       /* Enable fine timing estimation */
    bool foff_est_en;                     /* Enable frequency offset estimation/tracking */
    bool phase_est_en;                    /* Enable Phase offset estimation */

    complex float ** W;
    complex float * pilots;
    complex float * pilot_samples;
    complex float * rxbuf;
    float * w;
    kiss_fft_cfg sync_fft_cfg;      

    complex float * coarse_rxbuf;

    /* Demodulator data */

    complex float ** rx_sym;
    float * rx_amp;
    float * aphase_est_pilot_log;
};

#ifdef __cplusplus
}
#endif

#endif
