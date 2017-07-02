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

#ifndef M_PI
#define M_PI        3.14159265358979323846f  /* math constant */
#endif

#define TAU         (2.0f * M_PI)            /* mathematical constant */

#define OFDM_NC     16
#define OFDM_TS     0.018f
#define OFDM_RS     (1.0f / OFDM_TS)
#define OFDM_FS     8000.0f
#define OFDM_BPS    2
#define OFDM_TCP    0.002f
#define OFDM_NS     8
#define OFDM_CENTRE 1500.0f

/* To prevent C99 warning */

#define OFDM_M      144
#define OFDM_NCP    16

#ifdef OLD_STYLE
/* This will produce a warning in C99 as (int) makes these variable */

#define OFDM_M      ((int)(OFDM_FS / OFDM_RS))
#define OFDM_NCP    ((int)(OFDM_TCP * OFDM_FS))
#endif

#define OFDM_FTWINDOWWIDTH       11
#define OFDM_BITSPERFRAME        ((OFDM_NS - 1) * (OFDM_NC * OFDM_BPS))
#define OFDM_ROWSPERFRAME        (OFDM_BITSPERFRAME / (OFDM_NC * OFDM_BPS))
#define OFDM_SAMPLESPERFRAME     (OFDM_NS * (OFDM_M + OFDM_NCP))
#define OFDM_MAX_SAMPLESPERFRAME (OFDM_SAMPLESPERFRAME + (OFDM_M + OFDM_NCP)/4)
#define OFDM_RXBUF               (3 * OFDM_SAMPLESPERFRAME + 3 * (OFDM_M + OFDM_NCP))

struct OFDM {
    float foff_est_gain;
    float foff_est_hz;

    int verbose;
    int sample_point;
    int timing_est;
    int nin;

    bool timing_en;
    bool foff_est_en;
    bool phase_est_en;

    complex float pilot_samples[OFDM_M + OFDM_NCP];
    complex float W[OFDM_NC + 2][OFDM_M];
    complex float rxbuf[OFDM_RXBUF];
    complex float pilots[OFDM_NC + 2];
    float w[OFDM_NC + 2];
    
    /* Demodulator data */

    complex float rx_sym[OFDM_NS + 3][OFDM_NC + 2];
    complex float rx_np[OFDM_ROWSPERFRAME * OFDM_NC];
    float rx_amp[OFDM_ROWSPERFRAME * OFDM_NC];
    float aphase_est_pilot_log[OFDM_ROWSPERFRAME * OFDM_NC];
};

#ifdef __cplusplus
}
#endif

#endif
