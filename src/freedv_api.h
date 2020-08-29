/*---------------------------------------------------------------------------*\

  FILE........: freedv_api.h
  AUTHOR......: David Rowe
  DATE CREATED: August 2014

  Library of API functions that implement FreeDV "modes", useful for
  embedding FreeDV in other programs.  Please see:

  1. README_freedv.md
  2. Notes function use in freedv_api.c
  3. The sample freedv_tx.c and freedv_rx.c programs

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

#ifndef __FREEDV_API__
#define __FREEDV_API__

#if defined _WIN32 || defined __CYGWIN__
#ifdef __GNUC__
#define API_EXPORT __attribute__ ((dllexport))
#else
#define API_EXPORT __declspec(dllexport)
#endif
#else
#define API_EXPORT __attribute__ ((visibility ("default")))
#endif

// This declares a single-precision (float) complex number
#include <sys/types.h>

#include "comp.h"

#ifdef __cplusplus
  extern "C" {
#endif

// available modes
#define FREEDV_MODE_1600        0
#define FREEDV_MODE_2400A       3
#define FREEDV_MODE_2400B       4
#define FREEDV_MODE_800XA       5
#define FREEDV_MODE_700C        6
#define FREEDV_MODE_700D        7
#define FREEDV_MODE_2020        8

// Sample rates used
#define FREEDV_FS_8000          8000
#define FREEDV_FS_16000         16000

#ifndef FREEDV_MODE_EN_DEFAULT
#define FREEDV_MODE_EN_DEFAULT 1
#endif

// By default we enable all modes.  Disable during compile time e.g
// -DFREEDV_MODE_1600_EN=0 will enable all but FreeDV 1600.  Or the
// other way round -DFREEDV_MODE_EN_DEFAULT=0 -DFREEDV_MODE_1600_EN=1
// will enable only FreeDV 1600

#if !defined(FREEDV_MODE_1600_EN)
        #define FREEDV_MODE_1600_EN FREEDV_MODE_EN_DEFAULT
#endif
#if !defined(FREEDV_MODE_700C_EN)
        #define FREEDV_MODE_700C_EN FREEDV_MODE_EN_DEFAULT
#endif
#if !defined(FREEDV_MODE_700D_EN)
        #define FREEDV_MODE_700D_EN FREEDV_MODE_EN_DEFAULT
#endif
#if !defined(FREEDV_MODE_2400A_EN)
        #define FREEDV_MODE_2400A_EN FREEDV_MODE_EN_DEFAULT
#endif
#if !defined(FREEDV_MODE_2400B_EN)
        #define FREEDV_MODE_2400B_EN FREEDV_MODE_EN_DEFAULT
#endif
#if !defined(FREEDV_MODE_800XA_EN)
        #define FREEDV_MODE_800XA_EN FREEDV_MODE_EN_DEFAULT
#endif
#if !defined(FREEDV_MODE_2020_EN)
        #define FREEDV_MODE_2020_EN FREEDV_MODE_EN_DEFAULT
#endif

#define FDV_MODE_ACTIVE(mode_name, var)  ((mode_name##_EN) == 0 ? 0: (var) == mode_name)

// optional operator control of 700D state machine
#define FREEDV_SYNC_UNSYNC 0                 // force sync state machine to lose sync, and search for new sync
#define FREEDV_SYNC_AUTO   1                 // falls out of sync automatically
#define FREEDV_SYNC_MANUAL 2                 // fall out of sync only under operator control

// struct that hold state information for one freedv instance
struct freedv;

// Dummy structure for (currently) deprecated call
struct freedv_advanced {
    int interleave_frames;
};

// Called when text message char is decoded
typedef void (*freedv_callback_rx)(void *, char);
// Called when new text message char is needed
typedef char (*freedv_callback_tx)(void *);
typedef void (*freedv_calback_error_pattern)
       (void *error_pattern_callback_state, short error_pattern[], int sz_error_pattern);

// Protocol bits are packed MSB-first
// Called when a frame containing protocol data is decoded
typedef void (*freedv_callback_protorx)(void *, char *);
// Called when a frame containing protocol data is to be sent
typedef void (*freedv_callback_prototx)(void *, char *);

// Data packet callbacks
// Called when a packet has been received
typedef void (*freedv_callback_datarx)(void *, unsigned char *packet, size_t size);
// Called when a new packet can be send
typedef void (*freedv_callback_datatx)(void *, unsigned char *packet, size_t *size);


/*---------------------------------------------------------------------------*\
 
                                 FreeDV API functions

\*---------------------------------------------------------------------------*/

// open, close ----------------------------------------------------------------

API_EXPORT struct freedv *freedv_open_advanced(int mode, struct freedv_advanced *adv);
API_EXPORT struct freedv *freedv_open(int mode);
API_EXPORT void freedv_close   (struct freedv *freedv);

// Transmit -------------------------------------------------------------------

API_EXPORT void freedv_tx        (struct freedv *freedv, short mod_out[], short speech_in[]);
API_EXPORT void freedv_comptx    (struct freedv *freedv, COMP  mod_out[], short speech_in[]);
API_EXPORT void freedv_rawdatatx (struct freedv *f, short mod_out[], unsigned char *packed_payload_bits);
API_EXPORT void freedv_datatx    (struct freedv *f, short mod_out[]);
API_EXPORT int  freedv_data_ntxframes (struct freedv *freedv);

// Receive -------------------------------------------------------------------

API_EXPORT int freedv_nin       (struct freedv *freedv);
API_EXPORT int freedv_rx        (struct freedv *freedv, short speech_out[], short demod_in[]);
API_EXPORT int freedv_shortrx   (struct freedv *freedv, short speech_out[], short demod_in[], float gain);
API_EXPORT int freedv_floatrx   (struct freedv *freedv, short speech_out[], float demod_in[]);
API_EXPORT int freedv_comprx    (struct freedv *freedv, short speech_out[], COMP  demod_in[]);
API_EXPORT int freedv_rawdatarx (struct freedv *freedv, unsigned char *packed_payload_bits, short demod_in[]);

// Rawdata -------------------------------------------------------------------

API_EXPORT int freedv_codec_frames_from_rawdata(struct freedv *freedv, unsigned char *codec_frames, unsigned char *rawdata);
API_EXPORT int freedv_rawdata_from_codec_frames(struct freedv *freedv, unsigned char *rawdata, unsigned char *codec_frames);

// Set parameters ------------------------------------------------------------

API_EXPORT void freedv_set_callback_txt            (struct freedv *freedv, freedv_callback_rx rx, freedv_callback_tx tx, void *callback_state);
API_EXPORT void freedv_set_callback_protocol       (struct freedv *freedv, freedv_callback_protorx rx, freedv_callback_prototx tx, void *callback_state);
API_EXPORT void freedv_set_callback_data           (struct freedv *freedv, freedv_callback_datarx datarx, freedv_callback_datatx datatx, void *callback_state);
API_EXPORT void freedv_set_test_frames		(struct freedv *freedv, int test_frames);
API_EXPORT void freedv_set_test_frames_diversity	(struct freedv *freedv, int test_frames_diversity);
API_EXPORT void freedv_set_smooth_symbols		(struct freedv *freedv, int smooth_symbols);
API_EXPORT void freedv_set_squelch_en		(struct freedv *freedv, int squelch_en);
API_EXPORT void freedv_set_snr_squelch_thresh	(struct freedv *freedv, float snr_squelch_thresh);
API_EXPORT void freedv_set_clip	                (struct freedv *freedv, int val);
API_EXPORT void freedv_set_total_bit_errors    	(struct freedv *freedv, int val);
API_EXPORT void freedv_set_total_bits              (struct freedv *freedv, int val);
API_EXPORT void freedv_set_total_bit_errors_coded  (struct freedv *freedv, int val);
API_EXPORT void freedv_set_total_bits_coded        (struct freedv *freedv, int val);
API_EXPORT void freedv_set_callback_error_pattern  (struct freedv *freedv, freedv_calback_error_pattern cb, void *state);
API_EXPORT void freedv_set_varicode_code_num       (struct freedv *freedv, int val);
API_EXPORT void freedv_set_data_header             (struct freedv *freedv, unsigned char *header);
API_EXPORT void freedv_set_carrier_ampl            (struct freedv *freedv, int c, float ampl);
API_EXPORT void freedv_set_sync                    (struct freedv *freedv, int sync_cmd);
API_EXPORT void freedv_set_verbose                 (struct freedv *freedv, int verbosity);
API_EXPORT void freedv_set_tx_bpf                  (struct freedv *freedv, int val);
API_EXPORT void freedv_set_dpsk                    (struct freedv *freedv, int val);
API_EXPORT void freedv_set_ext_vco                 (struct freedv *f, int val);
API_EXPORT void freedv_set_phase_est_bandwidth_mode(struct freedv *f, int val);
API_EXPORT void freedv_set_eq                      (struct freedv *f, int val);
      
// Get parameters -------------------------------------------------------------------------

struct MODEM_STATS;

API_EXPORT int freedv_get_version(void);
API_EXPORT char *freedv_get_hash(void);
API_EXPORT int freedv_get_mode                 (struct freedv *freedv);
API_EXPORT void freedv_get_modem_stats         (struct freedv *freedv, int *sync, float *snr_est);
API_EXPORT void freedv_get_modem_extended_stats(struct freedv *freedv, struct MODEM_STATS *stats);
API_EXPORT int freedv_get_test_frames	    (struct freedv *freedv);

API_EXPORT int freedv_get_speech_sample_rate   (struct freedv *freedv);
API_EXPORT int freedv_get_n_speech_samples	    (struct freedv *freedv);
API_EXPORT int freedv_get_n_max_speech_samples (struct freedv *freedv);

API_EXPORT int freedv_get_modem_sample_rate    (struct freedv *freedv);
API_EXPORT int freedv_get_modem_symbol_rate    (struct freedv *freedv);
API_EXPORT int freedv_get_n_max_modem_samples  (struct freedv *freedv);
API_EXPORT int freedv_get_n_nom_modem_samples  (struct freedv *freedv);

// bit error rate stats
API_EXPORT int freedv_get_total_bits	    (struct freedv *freedv);
API_EXPORT int freedv_get_total_bit_errors	    (struct freedv *freedv);
API_EXPORT int freedv_get_total_bits_coded     (struct freedv *freedv);
API_EXPORT int freedv_get_total_bit_errors_coded(struct freedv *freedv);
API_EXPORT int freedv_get_uncorrected_errors   (struct freedv *freedv);

API_EXPORT int freedv_get_sync		    (struct freedv *freedv);
API_EXPORT int freedv_get_sync_interleaver	    (struct freedv *freedv);

// access to speech codec states
API_EXPORT struct FSK * freedv_get_fsk         (struct freedv *f);
API_EXPORT struct CODEC2 *freedv_get_codec2    (struct freedv *freedv);

API_EXPORT int freedv_get_bits_per_codec_frame (struct freedv *freedv);
API_EXPORT int freedv_get_bits_per_modem_frame (struct freedv *freedv);
API_EXPORT int freedv_get_sz_error_pattern     (struct freedv *freedv);
API_EXPORT int freedv_get_protocol_bits        (struct freedv *freedv);

#ifdef __cplusplus
}
#endif

#endif //__FREEDV_API__
