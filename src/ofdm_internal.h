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
#define M_PI        3.14159265358979323846f
#endif

#define TAU         (2.0f * M_PI)
#define ROT45       (M_PI / 4.0f)

#define cmplx(value) (cosf(value) + sinf(value) * I)
#define cmplxconj(value) (cosf(value) + sinf(value) * -I)

/* modem state machine states */
typedef enum {
    search,
    trial,
    synced
} State;

typedef enum {
    unsync,             /* force sync state machine to lose sync, and search for new sync */
    autosync,           /* falls out of sync automatically */
    manualsync          /* fall out of sync only under operator control */
} Sync;

/* phase estimator bandwidth options */

typedef enum {
    low_bw,             /* can only track a narrow freq offset, but accurate         */
    high_bw             /* can track wider freq offset, but less accurate at low SNR */
} PhaseEstBandwidth;

/*
 * User-defined configuration for OFDM modem.  Used to set up
 * constants at init time, e.g. for different bit rate modems.
 */

struct OFDM_CONFIG {
    float tx_centre; /* TX Centre Audio Frequency */
    float rx_centre; /* RX Centre Audio Frequency */
    float fs;  /* Sample Frequency */
    float rs;  /* Modulation Symbol Rate */
    float ts;  /* symbol duration */
    float tcp; /* Cyclic Prefix duration */
    float ofdm_timing_mx_thresh;

    int nc;  /* Number of carriers */
    int ns;  /* Number of Symbol frames */
    int bps;   /* Bits per Symbol */
    int txtbits; /* number of auxiliary data bits */
    int ftwindowwidth;
};

struct OFDM {
    // Pointers

    struct quisk_cfFilter *ofdm_tx_bpf;
    
    complex float *pilot_samples;
    complex float *rxbuf;
    complex float *pilots;
    complex float **rx_sym;
    complex float *rx_np;

    float *rx_amp;
    float *aphase_est_pilot_log;

    uint8_t *tx_uw;

    // State enums
    State sync_state;
    State last_sync_state;
    State sync_state_interleaver;
    State last_sync_state_interleaver;

    // Sync enums
    Sync sync_mode;

    // Phase enums
    PhaseEstBandwidth phase_est_bandwidth;

    int phase_est_bandwidth_mode;

    // Complex
    complex float foff_metric;
     
    // Float
    float foff_est_gain;
    float foff_est_hz;
    float timing_mx;
    float coarse_foff_est_hz;
    float timing_norm;
    float sig_var;
    float noise_var;
    float mean_amp;

    // Integer
    int clock_offset_counter;
    int verbose;
    int sample_point;
    int timing_est;
    int timing_valid;
    int nin;
    int uw_errors;
    int sync_counter;
    int frame_count;
    int frame_count_interleaver;

    // Boolean
    bool sync_start;
    bool sync_end;
    bool timing_en;
    bool foff_est_en;
    bool phase_est_en;
    bool tx_bpf_en;
    bool dpsk;
};

/* Prototypes */

complex float qpsk_mod(int *);
void qpsk_demod(complex float, int *);
void ofdm_txframe(struct OFDM *, complex float *, complex float []);
void ofdm_assemble_modem_frame(struct OFDM *, uint8_t [], uint8_t [], uint8_t []);
void ofdm_assemble_modem_frame_symbols(complex float [], COMP [], uint8_t []);
void ofdm_disassemble_modem_frame(struct OFDM *, uint8_t [], COMP [], float [], short []);
void ofdm_rand(uint16_t [], int);
void ofdm_generate_payload_data_bits(uint8_t [], int);
int ofdm_get_phase_est_bandwidth_mode(struct OFDM *);
void ofdm_set_phase_est_bandwidth_mode(struct OFDM *, int);

#ifdef __cplusplus
}
#endif

#endif

