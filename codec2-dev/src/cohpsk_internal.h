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
#define ND               2                           /* diversity factor ND 1 is no diveristy, ND we have orginal plus 
                                                        one copy */

#include "fdmdv_internal.h"
#include "kiss_fft.h"

struct COHPSK {
    float        pilot2[2*NPILOTSFRAME][COHPSK_NC];    
    float        phi_[NSYMROW][COHPSK_NC*ND];           /* phase estimates for this frame of rx data symbols     */
    float        amp_[NSYMROW][COHPSK_NC*ND];           /* amplitude estimates for this frame of rx data symbols */
    COMP         rx_symb[NSYMROW][COHPSK_NC*ND];        /* demodulated symbols                                   */
    kiss_fft_cfg fft_coarse_fest;
    float        f_est;
    COMP         rx_filter_memory[COHPSK_NC*ND][NFILTER];
    COMP         ct_symb_buf[NCT_SYMB_BUF][COHPSK_NC*ND];
    int          ct;                                 /* coarse timing offset in symbols                       */
    float        f_fine_est;
    COMP         ff_rect;
    COMP         ff_phase;
    COMP         ct_symb_ff_buf[NSYMROWPILOT+2][COHPSK_NC*ND];
    int          sync;
    int          sync_timer;

    struct FDMDV *fdmdv;
};


#endif
