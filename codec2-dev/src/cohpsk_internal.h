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

#define COARSE_FEST_NDFT 1024
#define NCT_SYMB_BUF     (2*NSYMROWPILOT+2)

#include "kiss_fft.h"

struct COHPSK {
    float        pilot2[2*NPILOTSFRAME][PILOTS_NC];    
    float        phi_[NSYMROW][PILOTS_NC];           /* phase estimates for this frame of rx data symbols     */
    float        amp_[NSYMROW][PILOTS_NC];           /* amplitude estimates for this frame of rx data symbols */
    COMP         rx_symb[NSYMROW][PILOTS_NC];        /* demodulated symbols                                   */
    kiss_fft_cfg fft_coarse_fest;
    float        f_est;
    COMP         ct_symb_buf[NCT_SYMB_BUF][PILOTS_NC];
    float        f_fine_est;
    int          ct;
    COMP         ff_rect;
};


#endif
