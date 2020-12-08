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

#ifndef __FREEDV_API_INTERNAL__
#define __FREEDV_API_INTERNAL__

#include "varicode.h"
#include "fsk.h"
#include "fmfsk.h"
#include "codec2_fdmdv.h"
#include "codec2_cohpsk.h"
#ifdef __LPCNET__
#include "lpcnet_freedv.h"
#endif
#include "freedv_api.h"

#ifdef __cplusplus
  extern "C" {
#endif

// Experimentally derived fudge factors to normalise Tx power across modes
#define NORM_PWR_COHPSK  1.74   
#define NORM_PWR_FSK     0.193 
#define NORM_PWR_OFDM    1.00

// identifiers for non Codec 2 Speech codecs, make sure no overlap with CODEC2_XXX modes
#define CODEC_MODE_LPCNET_1733 100

// Return code flags for freedv_*rx* functions
#define RX_TRIAL_SYNC       0x1       // demodulator has trial sync
#define RX_SYNC             0x2       // demodulator has sync
#define RX_BITS             0x4       // data bits have been returned
#define RX_BIT_ERRORS       0x8       // FEC may not have corrected all bit errors (not all parity checks OK)
                                      
extern char *rx_sync_flags_to_text[]; // converts flags above to more meaningful text
      
struct freedv {
    int                  mode;

    // states for various modules we support
    struct CODEC2       *codec2;
    struct FDMDV        *fdmdv;
    struct COHPSK       *cohpsk;
    struct FSK          *fsk;
    struct FMFSK        *fmfsk;
    struct OFDM         *ofdm;
    struct LDPC         *ldpc;
    struct MODEM_STATS   stats;
#ifdef __LPCNET__
    struct LPCNetFreeDV *lpcnet;
#endif
    
    struct freedv_vhf_deframer * deframer;      // Extracts frames from VHF stream

    struct quisk_cfFilter * ptFilter7500to8000; // Filters to change to/from 7500 and 8000 sps for 700 .... 700C
    struct quisk_cfFilter * ptFilter8000to7500;

    int                  n_speech_samples;       // number of speech samples we need for each freedv_tx() call
                                                 // num of speech samples output by freedv_rx() call
    int                  n_nom_modem_samples;    // size of tx modem sample buffers
    int                  n_max_modem_samples;    // make your rx modem sample buffers this big
    int                  n_nat_modem_samples;    // tx modem sample block length as used by the modem before interpolation to output
                                                 // usually the same as n_nom_modem_samples, except for 700C
    int                  modem_sample_rate;      // Caller is responsible for meeting this
    int                  modem_symbol_rate;      // Useful for ext_vco operation on 2400A and 800XA
    int                  speech_sample_rate;     // 8 kHz or 16 kHz (high fidelity)

    int                  bits_per_codec_frame;   // one of modem codec frames per modem frame
    int                  bits_per_modem_frame;   // number of modem payload bits in each modem frame (usually compressed speech)
    int                  n_codec_frames;         // number of codec frames in each modem frame
    uint8_t             *tx_payload_bits;        // payload bits (usually compressed speech) for a modem frame ...
    uint8_t             *rx_payload_bits;        // ... one bit per char for some modes, packed for others

    /* FDMDV buffers for FreeDV 1600 -------------------------------------------------------------*/
    
    int                 *fdmdv_bits;
    int                 *fdmdv_tx_bits;
    int                 *fdmdv_rx_bits;

    /* test frame states -------------------------------------------------------------------------*/
    
    int                 *ptest_bits_coh;
    int                 *ptest_bits_coh_end;

    int                  test_frames;            // set this baby for 1 to tx/rx test frames to look at bit error stats
    int                  test_frames_diversity;  // 1 -> used combined carriers for error counting on 700 waveforms
    int                  test_frame_sync_state;
    int                  test_frame_sync_state_upper;  // when test_frames_diveristy==0 we need extra states for upper carriers
    int                  test_frame_count;
    int                  total_bits;
    int                  total_bit_errors;
    int                  total_bits_coded;
    int                  total_bit_errors_coded;
    int                  sz_error_pattern;

    /* optional user defined function to pass error pattern when a test frame is received */

    void                *error_pattern_callback_state;
    void (*freedv_put_error_pattern)(void *error_pattern_callback_state, short error_pattern[], int sz_error_pattern);

    /* Misc ---------------------------------------------------------------------------------------------*/

    int                 *tx_bits;                            /* FSK modem frame under construction */
    int                  tx_sync_bit;
    int                  frames;
    int                  clip;                               /* non-zero for cohpsk modem output clipping for low PAPR */
    int                  sync;
    int                  evenframe;
    float                snr_est;
    float                snr_squelch_thresh;
    int                  squelch_en;
    int                  nin, nin_prev;
    int                  verbose;
    int                  ext_vco;                            /* 2400A/800XA use external VCO flag */
    float               *passthrough_2020;                   /* 2020 interpolating filter */
    float                tx_amp;                             /* amplitude of tx samples */
    
    int                  ofdm_bitsperframe;
    int                  ofdm_nuwbits;
    int                  ofdm_ntxtbits;
    int                  rx_status;
    
    /* Varicode txt channel states ----------------------------------------------------------------------*/
    
    struct VARICODE_DEC  varicode_dec_states;
    short                tx_varicode_bits[VARICODE_MAX_BITS];
    int                  nvaricode_bits;
    int                  varicode_bit_index;

    /* interleaved LDPC OFDM states ---------------------------------------------------------------------*/

    COMP                *codeword_symbols;
    float               *codeword_amps;
    int                  modem_frame_count_tx;       // modem frame counter for tx side
    int                  modem_frame_count_rx;       // modem frame counter for rx side
    COMP                *mod_out;                    // output buffer of intereaved frames
    
    /* user defined function ptrs to produce and consume ASCII
      characters using aux txt channel */

    char (*freedv_get_next_tx_char)(void *callback_state);
    void (*freedv_put_next_rx_char)(void *callback_state, char c);
    void                *callback_state;
    
    /* user defined functions to produce and consume protocol bits */
    /* Protocol bits are packed MSB-first */
    void (*freedv_put_next_proto)(void *callback_state, char *proto_bits_packed);
    void (*freedv_get_next_proto)(void *callback_state, char *proto_bits_packed);
    void *proto_callback_state;
    int n_protocol_bits;

    /* states needed for FSK LDPC */
    float   *frame_llr;
    int      frame_llr_size, frame_llr_nbits;
    float   *twoframes_llr;
    uint8_t *twoframes_hard;
    int      fsk_ldpc_thresh1, fsk_ldpc_thresh2, fsk_ldpc_baduw_thresh;
    int      fsk_ldpc_state,  fsk_ldpc_best_location, fsk_ldpc_baduw;
    float    fsk_ldpc_snr;
};

// open function for each mode
      
void freedv_1600_open(struct freedv *f);
void freedv_700c_open(struct freedv *f);
void freedv_700d_open(struct freedv *f);
void freedv_2020_open(struct freedv *f);
void freedv_2400a_open(struct freedv *f);
void freedv_2400b_open(struct freedv *f);
void freedv_800xa_open(struct freedv *f);
void freedv_fsk_ldpc_open(struct freedv *f, struct freedv_advanced *adv);

// each mode has tx and rx functions in various flavours for real and complex valued samples

void freedv_comptx_fdmdv_1600(struct freedv *f, COMP mod_out[]);
int freedv_comprx_fdmdv_1600(struct freedv *f, COMP demod_in[]);
      
void freedv_comptx_700c(struct freedv *f, COMP mod_out[]);
int freedv_comprx_700c(struct freedv *f, COMP demod_in_8kHz[]);

void freedv_comptx_700d(struct freedv *f, COMP mod_out[]);
int freedv_comp_short_rx_700d(struct freedv *f, void *demod_in_8kHz, int demod_in_is_short, float gain);

void freedv_comptx_2020(struct freedv *f, COMP mod_out[]);
int freedv_comprx_2020(struct freedv *f, COMP demod_in[]);

void freedv_comptx_fsk_voice(struct freedv *f, COMP mod_out[]);
void freedv_tx_fsk_voice(struct freedv *f, short mod_out[]);
void freedv_tx_fsk_data(struct freedv *f, short mod_out[]);
int freedv_comprx_fsk(struct freedv *f, COMP demod_in[]);
int freedv_floatrx(struct freedv *f, short speech_out[], float demod_in[]);

void freedv_tx_fsk_ldpc_data(struct freedv *f, COMP mod_out[]);
void freedv_tx_fsk_ldpc_data_preamble(struct freedv *f, COMP mod_out[], int npreamble_bits, int npreamble_samples);
int freedv_rx_fsk_ldpc_data(struct freedv *f, COMP demod_in[]);
      
int freedv_bits_to_speech(struct freedv *f, short speech_out[], short demod_in[], int rx_status);

#ifdef __cplusplus
}
#endif

#endif

