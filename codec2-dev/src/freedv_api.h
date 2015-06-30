/*---------------------------------------------------------------------------*\
                                                                             
  FILE........: freedv_api.h
  AUTHOR......: David Rowe
  DATE CREATED: August 2014
                                                                             
  Library of API functions that implement FreeDV "modes", useful for
  embedding FreeDV in other programs.  Please see the documentation
  for each function in freedv_api.c, and the sample freedv_tx.c and
  freedv_rx.c programs.
                                                                             
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

#define FREEDV_MODE_1600        0
#define FREEDV_MODE_700         1

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
    int                  n_nom_modem_samples;    // size of tx and most rx modenm sample buffers
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

    int                  sync;
    int                  evenframe;
    float                snr_est;
    float                snr_squelch_thresh;
    float                squelch_en;
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

struct freedv *freedv_open(int mode);
void freedv_close(struct freedv *freedv);

void freedv_tx(struct freedv *f, short mod_out[], short speech_in[]);
void freedv_comptx(struct freedv *f, COMP mod_out[], short speech_in[]);

int freedv_nin(struct freedv *f);
int freedv_rx(struct freedv *f, short speech_out[], short demod_in[]);
int freedv_floatrx(struct freedv *f, short speech_out[], float demod_in[]);
int freedv_comprx(struct freedv *f, short speech_out[], COMP demod_in[]);

#endif

#ifdef __cplusplus
}
#endif
