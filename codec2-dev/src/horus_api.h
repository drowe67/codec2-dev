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
      
#define HORUS_MODE_BINARY            0
#define HORUS_MODE_RTTY              1

struct horus;
struct MODEM_STATS;

struct horus *horus_open  (int mode);
void          horus_close (struct horus *hstates);

/* call before horus_rx() to determine how many shorts to pass in */

uint32_t      horus_nin   (struct horus *hstates);

/* returns 1 if ascii_out[] is valid */
      
int           horus_rx    (struct horus *hstates, char ascii_out[], short demod_in[]);

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

/* how much storage you need for demod_in[] and  ascii_out[] */
      
int           horus_get_max_demod_in         (struct horus *hstates);
int           horus_get_max_ascii_out_len    (struct horus *hstates);

#endif

#ifdef __cplusplus
}
#endif
