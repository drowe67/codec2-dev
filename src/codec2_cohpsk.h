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

#define COHPSK_BITS_PER_FRAME     56              /* hard coded for now */
#define COHPSK_NC                  7              /* hard coded for now */
#define COHPSK_SAMPLES_PER_FRAME 600
#define COHPSK_RS                 75
#define COHPSK_FS               7500              /* note this is a wierd value to get an integer oversampling rate */

#include "comp.h"
#include "codec2_fdmdv.h"

struct COHPSK;

struct COHPSK *cohpsk_create(void);
void cohpsk_destroy(struct COHPSK *coh);
void cohpsk_mod(struct COHPSK *cohpsk, COMP tx_fdm[], int tx_bits[]);
void cohpsk_demod(struct COHPSK *cohpsk, int rx_bits[], int *reliable_sync_bit, COMP rx_fdm[], int *nin_frame);

#endif
