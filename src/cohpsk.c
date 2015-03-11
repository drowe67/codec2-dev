/*---------------------------------------------------------------------------*\
                                                                             
  FILE........: cohpsk.c
  AUTHOR......: David Rowe
  DATE CREATED: March 2015
                                                                             
  Functions that implement a coherent PSK FDM modem.
        
  TODO:
    [ ] matching octave function bits_to_dqpsk_symbols()
    [ ] framework for test program and octave UT
    [ ] testframe used by both Octave and C
                                                               
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
#include "pilots_coh.h"

static COMP qpsk_mod[] = {
    { 1.0, 0.0},
    { 0.0, 1.0},
    { 0.0,-1.0},
    {-1.0, 0.0}
};
    
/*---------------------------------------------------------------------------*\
                                                                             
                               FUNCTIONS

\*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*\
                                                       
  FUNCTION....: bits_to_qpsk_symbols()	     
  AUTHOR......: David Rowe			      
  DATE CREATED: March 2015

  Maps bits to parallel DQPSK symbols and inserts pilot symbols.

\*---------------------------------------------------------------------------*/

void bits_to_qpsk_symbols(COMP tx_symb[][PILOTS_NC], int tx_bits[], int framesize)
{
    COMP tx_symb_data[NSYMROW][PILOTS_NC];
    int   i, r, c, p_r, data_r;
    short bits;

    assert(COHPSK_NC == PILOTS_NC);
    assert((NSYMROW*PILOTS_NC)/2 == framesize);
 
    /*
      Organise QPSK symbols into a NSYMBROWS rows by PILOTS_NC cols matrix,
      each column is a carrier, time flows down the rows......

      Note: the "& 0x1" prevents and non binary tx_bits[] screwing up
      our lives.  Call me defensive.
    */

    for(c=0; c<PILOTS_NC; c++) {
        for(r=0; r<NSYMROW; r++) {
            i = c*NSYMROW + r;
            bits = (tx_bits[2*i]&0x1)<<1 | (tx_bits[2*i+1]&0x1);          
            tx_symb_data[r][c] = qpsk_mod[bits];
        }
    }

    /* "push" in rows of Nc pilots, one row every NS data symbols */
            
    for (r=0, p_r=0, data_r=0; p_r<NPILOTSFRAME; p_r++) {

        /* row of pilots */

        for(c=0; c<PILOTS_NC; c++) {
            tx_symb[r][c].real = pilots_coh[p_r][c];
            tx_symb[r][c].imag = 0.0;
        }
        r++;

        /* NS rows of data symbols */

        for(i=0; i<NS; data_r++,r++) {
            for(c=0; c<PILOTS_NC; c++)
                tx_symb[r][c] = tx_symb_data[data_r][c];
        }
    }

    assert(p_r == NPILOTSFRAME);
    assert(data_r == NSYMROW);
    assert(r == NSYMROWPILOT);
}
