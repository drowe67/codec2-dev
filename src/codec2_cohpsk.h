/*---------------------------------------------------------------------------*\
                                                                             
  FILE........: cohpsk.h
  AUTHOR......: David Rowe
  DATE CREATED: March 2015
                                                                             
  Functions that implement a coherent PSK FDM modem.
                                                               
\*---------------------------------------------------------------------------*/

/*
  Copyright (C) 2015 David Rowe

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

#ifndef __CODEC2_COHPSK__
#define __CODEC2_COHPSK__

#define COHPSK_BITS_PER_FRAME  32             /* hard coded for now */
#define COHPSK_NC               4             /* hard coded for now */

#include "comp.h"
#include "codec2_fdmdv.h"

struct COHPSK;

struct COHPSK *cohpsk_create(void);
void cohpsk_destroy(struct COHPSK *coh);
void bits_to_qpsk_symbols(COMP tx_symb[][COHPSK_NC], int tx_bits[], int nbits);
void qpsk_symbols_to_bits(struct COHPSK *coh, int rx_bits[], COMP ct_symb_buf[][COHPSK_NC]);
void coarse_freq_offset_est(struct COHPSK *coh, struct FDMDV *fdmdv, COMP ch_fdm_frame[], int sync, int *next_sync);
void frame_sync_fine_timing_est(struct COHPSK *coh, COMP ch_symb[][COHPSK_NC], int sync, int *next_sync);
int sync_state_machine(sync, next_sync);

#endif
