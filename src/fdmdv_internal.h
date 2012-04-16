/*---------------------------------------------------------------------------*\
                                                                             
  FILE........: fdmdv_internal.h
  AUTHOR......: David Rowe                                                          
  DATE CREATED: April 16 2012
                                                                             
  Header file for FDMDV internal functions, exposed via this header
  file for testing.
                                                                             
\*---------------------------------------------------------------------------*/

/*
  Copyright (C) 2012 David Rowe

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

#ifndef __FDMDV_INTERNAL__
#define __FDMDV_INTERNAL__

#include "comp.h"

void bits_to_dqpsk_symbols(COMP tx_symbols[], COMP prev_tx_symbols[], int tx_bits[], int *pilot_bit);

#endif
