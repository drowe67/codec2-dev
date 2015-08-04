/*---------------------------------------------------------------------------*\
                                                                             
  FILE........: freedv_api_internal.h
  AUTHOR......: David Rowe
  DATE CREATED: August 2014
                                                                             
  This declares the structure freedv.  A pointer to this structure is
  returned by the FreeDV API freedv_open() function.  The pointer is used
  by the other FreeDV API functions declared in freedv_api.h.  This
  structure is intended to be internal to the FreeDV API.  The public
  functions are declared in freedv_api.h.  Changes to this structure
  are expected.  Changes (except additions) to freedv_api.h are
  discouraged.
                                                                             
\*---------------------------------------------------------------------------*/

/*
  Copyright (C) 2014 David Rowe

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

#ifdef __cplusplus
  extern "C" {
#endif

#ifndef __FREEDV__

#include "varicode.h"
#include "codec2_fdmdv.h"
#include "codec2_cohpsk.h"

struct freedv {
    int                  mode;

    struct CODEC2       *codec2;
    struct FDMDV        *fdmdv;
    struct MODEM_STATS   stats;
    struct COHPSK       *cohpsk;

    int                  n_speech_samples;
    int                  n_nom_modem_samples;    // size of tx and most rx modem sample buffers
    int                  n_max_modem_samples;    // make your rx modem sample buffers this big

    int                  modem_sample_rate;      // ATM caller is responsible for meeting this (TBC)
    int                  clip;                   // non-zero for cohpsk modem output clipping for low PAPR

    unsigned char       *packed_codec_bits;
    int                 *codec_bits;
    int                 *tx_bits;
    int                 *fdmdv_bits;
    int                 *rx_bits;
    int                  tx_sync_bit;
    int                  smooth_symbols;
    float               *prev_rx_bits;

    int                 *ptest_bits_coh;
    int                 *ptest_bits_coh_end;

    int                  test_frames;            // set this baby for 1 to tx/rx test frames to look at bit error stats
    int                  test_frame_sync_state;
    int                  test_frame_count;
    int                  total_bits;
    int                  total_bit_errors;
    int                  sz_error_pattern;

    /* optional user defined function to pass error pattern when a test frame is received */

    void                *error_pattern_callback_state;
    void (*freedv_put_error_pattern)(void *error_pattern_callback_state, short error_pattern[], int sz_error_pattern);

    int                  sync;
    int                  evenframe;
    float                snr_est;
    float                snr_squelch_thresh;
    int                  squelch_en;
    int                  nin;

    struct VARICODE_DEC  varicode_dec_states;
    short                tx_varicode_bits[VARICODE_MAX_BITS];
    int                  nvaricode_bits;
    int                  varicode_bit_index;
    
    /* user defined function ptrs to produce and consume ASCII
      characters using aux txt channel */

    char (*freedv_get_next_tx_char)(void *callback_state);
    void (*freedv_put_next_rx_char)(void *callback_state, char c);

    void                *callback_state;

};

#endif

#ifdef __cplusplus
}
#endif
