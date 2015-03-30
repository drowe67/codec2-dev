/*---------------------------------------------------------------------------*\
                                                                             
  FILE........: cohpsk.c
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

/*---------------------------------------------------------------------------*\
                                                                             
                               INCLUDES

\*---------------------------------------------------------------------------*/

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "codec2_cohpsk.h"
#include "test_bits.h"
#include "cohpsk_defs.h"
#include "cohpsk_internal.h"
#include "pilots_coh.h"
#include "comp_prim.h"

static COMP qpsk_mod[] = {
    { 1.0, 0.0},
    { 0.0, 1.0},
    { 0.0,-1.0},
    {-1.0, 0.0}
};
    
/*---------------------------------------------------------------------------*\
                                                                             
                               FUNCTIONS

\*---------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------* \
                                                       
  FUNCTION....: cohpsk_create	     
  AUTHOR......: David Rowe			      
  DATE CREATED: Marcg 2015 

  Create and initialise an instance of the modem.  Returns a pointer
  to the modem states or NULL on failure.  One set of states is
  sufficient for a full duplex modem.

\*---------------------------------------------------------------------------*/

struct COHPSK *cohpsk_create(void)
{
    struct COHPSK *coh;
    int            r,c,p;

    coh = (struct COHPSK*)malloc(sizeof(struct COHPSK));
    if (coh == NULL)
        return NULL;

    /* set up buffer of tx pilot symbols for rx */

    for(r=0; r<2*NPILOTSFRAME; ) {
        for(p=0; p<NPILOTSFRAME; r++, p++) {
            for(c=0; c<PILOTS_NC; c++) {
                coh->pilot2[r][c] = pilots_coh[p][c];
            }
        }
    }

    return coh;
}


/*---------------------------------------------------------------------------*\
                                                       
  FUNCTION....: cohpsk_destroy	     
  AUTHOR......: David Rowe			      
  DATE CREATED: March 2015

  Destroy an instance of the modem.

\*---------------------------------------------------------------------------*/

void cohpsk_destroy(struct COHPSK *coh)
{
    assert(coh != NULL);
    free(coh);
}


/*---------------------------------------------------------------------------*\
                                                       
  FUNCTION....: bits_to_qpsk_symbols()	     
  AUTHOR......: David Rowe			      
  DATE CREATED: March 2015

  Rate Rs modulator.  Maps bits to parallel DQPSK symbols and inserts pilot symbols.

\*---------------------------------------------------------------------------*/

void bits_to_qpsk_symbols(COMP tx_symb[][PILOTS_NC], int tx_bits[], int nbits)
{
    int   i, r, c, p_r, data_r;
    short bits;

    assert(COHPSK_NC == PILOTS_NC);
    assert((NSYMROW*PILOTS_NC)*2 == nbits);
 
    /*
      Insert two rows of Nc pilots at beginning of data frame.

      Organise QPSK symbols into a NSYMBROWS rows by PILOTS_NC cols matrix,
      each column is a carrier, time flows down the cols......

      Note: the "& 0x1" prevents and non binary tx_bits[] screwing up
      our lives.  Call me defensive.
    */

    r = 0;
    for(p_r=0; p_r<2; p_r++) {
        for(c=0; c<PILOTS_NC; c++) {
            tx_symb[r][c].real = pilots_coh[p_r][c];
            tx_symb[r][c].imag = 0.0;
        }
        r++;
    }
    for(data_r=0; data_r<NSYMROW; data_r++, r++) {
        for(c=0; c<PILOTS_NC; c++) {
            i = c*NSYMROW + data_r;
            bits = (tx_bits[2*i]&0x1)<<1 | (tx_bits[2*i+1]&0x1);          
            tx_symb[r][c] = qpsk_mod[bits];
        }
    }
    
    assert(p_r == NPILOTSFRAME);
    assert(r == NSYMROWPILOT);
}


/*---------------------------------------------------------------------------*\
                                                       
  FUNCTION....: qpsk_symbols_to_bits()	     
  AUTHOR......: David Rowe			      
  DATE CREATED: March 2015

  Rate Rs demodulator. Extract pilot symbols and estimate amplitude and phase
  of each carrier.  Correct phase of data symbols and convert to bits.

  Further improvement.  In channels with rapidly changing phase by
  moderate Eb/No, we could perhaps do better by interpolating the
  phase across symbols rather than using the same phi_ for all symbols.

\*---------------------------------------------------------------------------*/

void qpsk_symbols_to_bits(struct COHPSK *coh, int rx_bits[], COMP ct_symb_buf[][COHPSK_NC])
{
    int   r, c, i;
    COMP  corr, rot, pi_on_4;
    float mag, phi_, amp_;
    short sampling_points[] = {1, 2, 7, 8};

    pi_on_4.real = cosf(M_PI/4); pi_on_4.imag = sinf(M_PI/4);
   
    /* Average pilots to get phase and amplitude estimates we assume
       there are two pilots at the start of each frame and two at the
       end */

    for(c=0; c<PILOTS_NC; c++) {
        corr.real = 0.0; corr.imag = 0.0; mag = 0.0;
        for(r=0; r<2*NPILOTSFRAME; r++) {
            corr = cadd(corr, fcmult(coh->pilot2[r][c], ct_symb_buf[sampling_points[r]][c]));
            mag  += cabsolute(ct_symb_buf[sampling_points[r]][c]);
        }
      
        phi_ = atan2f(corr.imag, corr.real);
        amp_ =  mag/2*NPILOTSFRAME;
        for(r=0; r<2*NPILOTSFRAME; r++) {
            coh->phi_[r][c] = phi_;
            coh->amp_[r][c] = amp_;
        }
    }

    /* now correct phase of data symbols and make decn on bits */

    for(c=0; c<PILOTS_NC; c++) {
        rot.real = cosf(coh->phi_[0][c]); rot.imag = -sinf(coh->phi_[0][c]);
        for (r=0; r<NSYMROW; r++) {
            i = c*NSYMROW + r;
            coh->rx_symb[r][c] = cmult(ct_symb_buf[NPILOTSFRAME + r][c], rot);
            rot = cmult(coh->rx_symb[r][c], pi_on_4);
            rx_bits[2*i+1] = rot.real < 0;
            rx_bits[2*i]   = rot.imag < 0;
        }
    }

}
