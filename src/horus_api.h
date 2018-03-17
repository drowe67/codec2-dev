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

#define HORUS_MODE_RTTY        0
#define HORUS_MODE_BINARY      1

struct MODEM_STATS;

struct horus *horus_open  (int mode);
void          horus_close (struct horus *hstates);

int           horus_nin   (struct horus *hstates);
int           horus_rx    (struct horus *hstates, char telemetry_out[], short demod_in[]);

int           horus_get_version              (void);
int           horus_get_mode                 (struct horus *hstates);
void          horus_get_modem_stats          (struct horus *hstates, int *sync, float *snr_est);
void          horus_get_modem_extended_stats (struct horus *hstates, struct MODEM_STATS *stats);

#endif

#ifdef __cplusplus
}
#endif
