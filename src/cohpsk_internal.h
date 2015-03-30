/*---------------------------------------------------------------------------*\
                                                                             
  FILE........: cohpsk_internal.h
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

#ifndef __COHPSK_INTERNAL__
#define __COHPSK_INTERNAL__

struct COHPSK {
    float pilot2[2*NPILOTSFRAME][PILOTS_NC];    
    float phi_[NSYMROW][PILOTS_NC];                     /* phase estimates for this frame of rx data symbols     */
    float amp_[NSYMROW][PILOTS_NC];                     /* amplitude estimates for this frame of rx data symbols */
    COMP  rx_symb[NSYMROW][PILOTS_NC];                  /* demodulated symbols                                   */
};

#endif
