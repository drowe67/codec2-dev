/*---------------------------------------------------------------------------*\

  FILE........: horus_api.h
  AUTHOR......: David Rowe
  DATE CREATED: March 2018

  Library of API functions that implement High Altitude Balloon (HAB)
  telemetry modems and protocols.

\*---------------------------------------------------------------------------*/

/*
  Copyright (C) 2018 David Rowe

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

#ifndef __HORUS_API__

#include <stdint.h>
#include "modem_stats.h"
      
/* Horus API Modes */
#define HORUS_MODE_BINARY            0  // Legacy binary mode
#define HORUS_MODE_BINARY_256BIT     1  // New 256-bit LDPC-encoded mode
#define HORUS_MODE_BINARY_128BIT     2  // New 128-bit LDPC-encoded mode
#define HORUS_MODE_RTTY              99 // RTTY Decoding


// Settings for Legacy Horus Binary Mode (Golay Encoding)
#define HORUS_BINARY_NUM_BITS                   360   /* fixed number of coded bits */
#define HORUS_BINARY_NUM_PAYLOAD_BYTES          22    /* fixed number of bytes in binary payload     */
#define HORUS_BINARY_DEFAULT_BAUD               100   /* Default baud rate of 100 baud */
#define HORUS_BINARY_DEFAULT_TONE_SPACING       270   /* Default tone spacing of 270 Hz */

// Settings for Horus Binary 256-bit mode (LDPC Encoding, r=1/3)
#define HORUS_BINARY_256BIT_NUM_BITS            768
#define HORUS_BINARY_256BIT_NUM_PAYLOAD_BYTES   32
#define HORUS_BINARY_256BIT_DEFAULT_BAUD        100   /* Default baud rate of 100 baud */
#define HORUS_BINARY_256BIT_DEFAULT_TONE_SPACING 270   /* Default tone spacing of 270 Hz */

// Settings for Horus Binary 128-bit mode (LDPC Encoding, r=1/3)
#define HORUS_BINARY_128BIT_NUM_BITS            384
#define HORUS_BINARY_128BIT_NUM_PAYLOAD_BYTES   16
#define HORUS_BINARY_128BIT_DEFAULT_BAUD        50   /* Default baud rate of 100 baud */
#define HORUS_BINARY_128BIT_DEFAULT_TONE_SPACING 270   /* Default tone spacing of 270 Hz */

// Settings for RTTY Decoders
#define HORUS_RTTY_NUM_BITS                     1000 /* Maximum 1000 bits in a RTTY frame (100 characters) */
#define HORUS_RTTY_DEFAULT_BAUD                 100 /* Default baud rate of 100 baud */

struct horus;
struct MODEM_STATS;

/*
 * Create an Horus Demod config/state struct using default mode parameters.
 * 
 * int mode - Horus Mode Type (refer list above)
 */
struct horus *horus_open  (int mode);

/*
 * Create an Horus Demod config/state struct with more customizations.
 * 
 * int mode - Horus Mode Type (refer list above)
 * int Rs - Symbol Rate (Hz)
 * int tx_tone_spacing - FSK Tone Spacing, to configure mask estimator. Set to -1 to disable mask estimator.
 */

struct horus *horus_open_advanced (int mode, int Rs, int tx_tone_spacing);

/*
 * Close a Horus demodulator struct and free memory.
 */
void          horus_close (struct horus *hstates);

/* call before horus_rx() to determine how many shorts to pass in */

uint32_t      horus_nin   (struct horus *hstates);

/*
 * Demodulate some number of Horus modem samples. The number of samples to be 
 * demodulated can be found by calling horus_nin().
 * 
 * Returns 1 if the data in ascii_out[] is valid.
 * 
 * struct horus *hstates - Horus API config/state struct, set up by horus_open / horus_open_advanced
 * char ascii_out[] - Buffer for returned packet / text.
 * short fsk_in[] - nin samples of modulated FSK.
 * int quadrature - Set to 1 if input samples are complex samples.
 */
      
int           horus_rx    (struct horus *hstates, char ascii_out[], short demod_in[], int quadrature);

/* set verbose level */
      
void horus_set_verbose(struct horus *hstates, int verbose);
      
/* functions to get information from API  */
      
int           horus_get_version              (void);
int           horus_get_mode                 (struct horus *hstates);
int           horus_get_Fs                   (struct horus *hstates);      
int           horus_get_mFSK                 (struct horus *hstates);      
void          horus_get_modem_stats          (struct horus *hstates, int *sync, float *snr_est);
void          horus_get_modem_extended_stats (struct horus *hstates, struct MODEM_STATS *stats);
int           horus_crc_ok                   (struct horus *hstates);
int           horus_get_total_payload_bits   (struct horus *hstates);
void          horus_set_total_payload_bits   (struct horus *hstates, int val);
void          horus_set_freq_est_limits      (struct horus *hstates, float fsk_lower, float fsk_upper);

/* how much storage you need for demod_in[] and  ascii_out[] */
      
int           horus_get_max_demod_in         (struct horus *hstates);
int           horus_get_max_ascii_out_len    (struct horus *hstates);

#endif

#ifdef __cplusplus
}
#endif
