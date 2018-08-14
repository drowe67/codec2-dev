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
  it under the terms of the GNU Lesser General Public License version 2.1, as
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
extern "C"
{
#endif

#include <complex.h>
#include <stdbool.h>
#include <stdint.h>

#include "codec2_ofdm.h"
#include "filter.h"

#ifndef M_PI
#define M_PI        3.14159265358979323846f  /* math constant */
#endif

#define TAU         (2.0f * M_PI)            /* mathematical constant */

/*
 * Contains user configuration for OFDM modem
 */

struct OFDM_CONFIG {
    float centre; /* Centre Audio Frequency */
    float fs;  /* Sample Frequency */
    float rs;  /* Baud Rate */
    float ts;  /* symbol duration */
    float tcp; /* Cyclic Prefix duration */
    float ofdm_timing_mx_thresh;

    int nc;  /* Number of carriers */
    int ns;  /* Number of Symbol frames */
    int bps;   /* Bits per Symbol */
    int txtbits; /* number of auxiliary data bits */
    int state_str; /* state string length */
    int ftwindowwidth;
};

struct OFDM {
    complex float *pilot_samples;
    complex float *rxbuf;
    complex float *pilots;
    complex float **rx_sym;
    complex float *rx_np;

    float *w;
    float *rx_amp;
    float *aphase_est_pilot_log;

    int *tx_uw;
    
    char *sync_state;
    char *last_sync_state;
    char *sync_state_interleaver;
    char *last_sync_state_interleaver;

    struct quisk_cfFilter ofdm_tx_bpf;
    
    complex float foff_metric;
    
    float foff_est_gain;
    float foff_est_hz;
    float timing_mx;
    float coarse_foff_est_hz;
    float timing_norm;
    float sig_var;
    float noise_var;
    float mean_amp;

    int clock_offset_counter;
    int verbose;
    int sample_point;
    int timing_est;
    int timing_valid;
    int nin;
    int uw_errors;
    int sync_counter;
    int frame_count;
    int sync_start;
    int sync_end;
    int sync_mode;
    int frame_count_interleaver;
    
    bool timing_en;
    bool foff_est_en;
    bool phase_est_en;
    bool tx_bpf_en;
};

/* function headers exposed for LDPC work */

complex float qpsk_mod(int *);
void qpsk_demod(complex float, int *);
void ofdm_txframe(struct OFDM *, complex float *, complex float []);
void ofdm_assemble_modem_frame(uint8_t [], uint8_t [], uint8_t []);
void ofdm_assemble_modem_frame_symbols(complex float [], COMP [], uint8_t []);
void ofdm_disassemble_modem_frame(struct OFDM *, int [], COMP [], float [], short []);
void ofdm_rand(uint16_t [], int);
void ofdm_generate_payload_data_bits(int payload_data_bits[], int data_bits_per_frame);

#ifdef __cplusplus
}
#endif

#endif
