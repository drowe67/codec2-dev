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

#define COHPSK_BITS_PER_FRAME 160             /* hard coded for now */
#define COHPSK_NC               4             /* hard coded for now */

#include "comp.h"

void bits_to_qpsk_symbols(COMP tx_symb[][COHPSK_NC], int tx_bits[], int framesize);

#endif
