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

#include <stdbool.h>
#include <stdint.h>

#include "compiler.h"
#include "codec2_ofdm.h"
#include "filter.h"

#ifdef __cplusplus
extern "C"
{
#endif

#ifndef M_PI
#define M_PI        3.14159265358979323846f
#endif

#define TAU         (2.0f * M_PI)
#define ROT45       (M_PI / 4.0f)

#define cmplx(value) complexf(cosf(value), sinf(value))
#define cmplxconj(value) complexf(cosf(value), -sinf(value))

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
    float tx_centre;   /* TX Centre Audio Frequency */
    float rx_centre;   /* RX Centre Audio Frequency */
    float fs;          /* Sample Frequency */
    float rs;          /* Symbol Rate */
    float ts;          /* symbol duration */
    float tcp;         /* Cyclic Prefix duration */
    float timing_mx_thresh;

    int nc;            /* Number of carriers */
    int ns;            /* Number of Symbol frames */
    int np;            /* number of modem frames per packet */
    int bps;           /* Bits per Symbol */
    int txtbits;       /* number of auxiliary data bits */
    int nuwbits;       /* number of unique word bits */
    int bad_uw_errors;
    int ftwindowwidth;
    int data_mode;     /* non-zero if this is a data mode */
    char *codename;    /* name of LDPC code used with this mode */
};

struct OFDM {
    struct OFDM_CONFIG config;
    
    /*
     * See 700D Part 4 Acquisition blog post and ofdm_dev.m routines
     * for how this was set
     */
    float timing_mx_thresh;
    
    int nc;
    int ns;	/* NS-1 = data symbols between pilots  */
    int bps; 	/* Bits per symbol */
    int m; 	/* duration of each symbol in samples */
    int ncp; 	/* duration of CP in samples */
    int np;     /* number of modem frames per packet. In some modes we want */
                /* the total packet of data to span multiple modem frames, e.g. HF data */
                /* and/or when the FEC codeword is larger than the one */
                /* modem frame.  In other modes (e.g. 700D/2020) Np=1, ie the modem frame */
                /* is the same length as the packet/FEC frame. */
    int ftwindowwidth;
    int bitsperframe;      /* total bits in all data symbols in modem frame */
    int bitsperpacket;     /* total bits in all data symbols in a packet */
    int rowsperframe;
    int samplespersymbol;
    int samplesperframe;
    int max_samplesperframe;
    int nuwframes;
    int nrxbuf;
    int ntxtbits;         /* reserve bits/frame for aux text information */
    int nuwbits;          /* number of unique word bits used to achieve packet frame sync */
    int bad_uw_errors;

    float tx_centre; /* TX Center frequency */
    float rx_centre; /* RX Center frequency */
    float fs; /* Sample rate */
    float ts; /* Symbol cycle time */
    float rs; /* Symbol rate */
    float tcp; /* Cyclic prefix duration */
    float tpacket; /* time for one packet in ms */
    float inv_m; /* 1/m */
    float tx_nlower; /* TX lowest carrier freq */
    float rx_nlower; /* RX lowest carrier freq */
    float doc; /* division of radian circle */
    
    // Pointers

    struct quisk_cfFilter *tx_bpf;
    
    complexf_t *pilot_samples;
    complexf_t *rxbuf;
    complexf_t *pilots;
    complexf_t **rx_sym;
    complexf_t *rx_np;
    complexf_t *tx_uw_syms;
    
    float *rx_amp;
    float *aphase_est_pilot_log;

    uint8_t *tx_uw;
    int *uw_ind;
    int *uw_ind_sym;

    // State enums
    State sync_state;
    State last_sync_state;

    // Sync enums
    Sync sync_mode;

    // Phase enums
    PhaseEstBandwidth phase_est_bandwidth;

    int phase_est_bandwidth_mode;

    // Complex
    complexf_t foff_metric;
     
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
    int modem_frame; /* increments for every modem frame in packet */
    int data_mode;
    
    // Boolean
    bool sync_start;
    bool sync_end;
    bool timing_en;
    bool foff_est_en;
    bool phase_est_en;
    bool tx_bpf_en;
    bool dpsk_en;

    char *codename;
};

/* Prototypes */

complexf_t qpsk_mod(int *);
complexf_t qam16_mod(int *);
void qpsk_demod(complexf_t, int *);
void qam16_demod(complexf_t, int *);
void ofdm_txframe(struct OFDM *, complexf_t *, complexf_t []);
void ofdm_assemble_qpsk_modem_packet(struct OFDM *, uint8_t [], uint8_t [], uint8_t []);
void ofdm_assemble_qpsk_modem_packet_symbols(struct OFDM *, complexf_t [], COMP [], uint8_t []);
void ofdm_disassemble_qpsk_modem_packet(struct OFDM *, complexf_t rx_syms[], float rx_amps[], COMP [], float [], short []);
void ofdm_extract_uw(struct OFDM *ofdm, complexf_t rx_syms[], float rx_amps[], uint8_t rx_uw[]);
void ofdm_rand(uint16_t [], int);
void ofdm_generate_payload_data_bits(uint8_t [], int);
int ofdm_get_phase_est_bandwidth_mode(struct OFDM *);
void ofdm_set_phase_est_bandwidth_mode(struct OFDM *, int);

#ifdef __cplusplus
}
#endif

#endif
