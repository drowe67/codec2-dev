/*---------------------------------------------------------------------------*\

  FILE........: codec2_ofdm.h
  AUTHORS.....: David Rowe & Steve Sampson
  DATE CREATED: June 2017

  External user references to the modem library.

\*---------------------------------------------------------------------------*/

/*
  Copyright (C) 2017 David Rowe

  All rights reserved.

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU Lesser General Public License version 2, as
  published by the Free Software Foundation.  This program is
  distributed in the hope that it will be useful, but WITHOUT ANY
  WARRANTY; without even the implied warranty of MERCHANTABILITY or
  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public
  License for more details.

  You should have received a copy of the GNU Lesser General Public License
  along with this program; if not, see <http://www.gnu.org/licenses/>.
*/

#ifndef CODEC2_OFDM_H
#define CODEC2_OFDM_H

#ifdef __cplusplus
extern "C" {
#endif

/* Includes */
    
#include <complex.h>
#include <stdbool.h>
    
#include "comp.h"

/* Defines */

#define OFDM_AMP_SCALE (2E5*1.1491/1.06)   /* use to scale to 16 bit short */
#define OFDM_CLIP (32767*0.35)             /* experimentally derived constant to reduce PAPR to about 8dB */
    
struct OFDM;

/* Default configuration for '700D' mode */
const struct OFDM_CONFIG * OFDM_CONFIG_700D;

/* create and destroy modem states */

struct OFDM *ofdm_create(const struct OFDM_CONFIG * config);
void ofdm_destroy(struct OFDM *);

/* signal processing */

void ofdm_mod(struct OFDM *, COMP *, const int *);
void ofdm_demod(struct OFDM *, int *, COMP *);
int  ofdm_sync_search(struct OFDM *ofdm, COMP *rxbuf_in);
void ofdm_sync_state_machine(struct OFDM *ofdm, int *rx_uw);

/* getters */
    
int ofdm_get_nin(struct OFDM *);
int ofdm_get_samples_per_frame(void);
int ofdm_get_max_samples_per_frame(void);
int ofdm_get_bits_per_frame(struct OFDM *);

/* option setters */

void ofdm_set_verbose(struct OFDM *, int);
void ofdm_set_timing_enable(struct OFDM *, bool);
void ofdm_set_foff_est_enable(struct OFDM *, bool);
void ofdm_set_phase_est_enable(struct OFDM *, bool);
void ofdm_set_off_est_hz(struct OFDM *, float);

#ifdef __cplusplus
}
#endif

#endif

