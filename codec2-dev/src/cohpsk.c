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

    coh = (struct COHPSK*)malloc(sizeof(struct COHPSK));
    if (coh == NULL)
        return NULL;

    /* set up buffer of 3 frames of tx pilot symbols */

    memcpy(&coh->tx_pilot_buf[0][0], pilots_coh, sizeof(pilots_coh));
    memcpy(&coh->tx_pilot_buf[NPILOTSFRAME][0], pilots_coh, sizeof(pilots_coh));
    memcpy(&coh->tx_pilot_buf[2*NPILOTSFRAME][0], pilots_coh, sizeof(pilots_coh));

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
    COMP tx_symb_data[NSYMROW][PILOTS_NC];
    int   i, r, c, p_r, data_r;
    short bits;

    assert(COHPSK_NC == PILOTS_NC);
    assert((NSYMROW*PILOTS_NC)*2 == nbits);
 
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

        for(i=0; i<NS; data_r++,r++,i++) {
            for(c=0; c<PILOTS_NC; c++) {
                tx_symb[r][c] = tx_symb_data[data_r][c];
            }
        }
    }

    assert(p_r == NPILOTSFRAME);
    assert(data_r == NSYMROW);
    assert(r == NSYMROWPILOT);
}


/*---------------------------------------------------------------------------*\
                                                       
  FUNCTION....: bits_to_qpsk_symbols()	     
  AUTHOR......: David Rowe			      
  DATE CREATED: March 2015

  Rate Rs demodulator. Extract pilot symbols and estimate amplitude and phase
  of each carrier.  Correct phase of data symbols and convert to bits.

\*---------------------------------------------------------------------------*/

void qpsk_symbols_to_bits(struct COHPSK *coh, int rx_bits[], COMP rx_symb[][COHPSK_NC])
{
    int   st, en, r, c, i, p_r, data_r;
    COMP  ch_est, rot, pi_on_4;

    pi_on_4.real = cosf(M_PI/4); pi_on_4.imag = sinf(M_PI/4);

    /* extract pilot and data symbols into 3 frame buffers */

    for (r=0, p_r=2*NPILOTSFRAME, data_r=2*NSYMROW; r<NSYMROWPILOT; p_r++) {

        printf("r: %d data_r: %d\n", r, data_r);
        /* copy row of pilots onto end of pilot buffer */

        for(c=0; c<PILOTS_NC; c++) {
            coh->rx_pilot_buf[p_r][c] = rx_symb[r][c];
            //printf("  %d %d %f %f", p_r, c, coh->rx_pilot_buf[p_r][c].real, coh->rx_pilot_buf[p_r][c].imag);
        }
        //printf("\n");
        r++;

        //printf("r: %d data_r: %d\n", r, data_r);
        /* copy NS rows of data symbols onto end of data symbol buffer */

        for(i=0; i<NS; data_r++,r++,i++) {
            for(c=0; c<PILOTS_NC; c++)
                coh->rx_symb_buf[data_r][c] = rx_symb[r][c];
        }
    }

    /* estimate channel amplitude and phase and correct data symbols in middle of buffer */

    for (r=0, data_r=NSYMROW; r<NSYMROW; r++, data_r) {

        /* pilots to use for correcting data_r-th symbol */

        st = NPILOTSFRAME + floor(r/NS) - floor(NP/2) + 1;
        en = st + NP - 1;
        assert(st >= 0);
        assert(en < 3*NPILOTSFRAME);

        printf("r: %d data_r: %d st: %d en: %d\n", r, data_r, st, en);

        /* iterate over all of the carriers */

        for (c=0; c<PILOTS_NC; c++) {

            /* estimate channel */

            ch_est.real = 0.0; ch_est.imag = 0.0;
            for (i=st; i<=en; i++)
                ch_est = cadd(ch_est, fcmult(coh->tx_pilot_buf[i][c], coh->rx_pilot_buf[i][c]));
            ch_est = fcmult(1.0/NP, ch_est);
            coh->phi_[r][c] = atan2(ch_est.imag,ch_est.real);
            coh->amp_[r][c] = cabsolute(ch_est);

            /* correct phase */

            rot.real = cosf(coh->phi_[r][c]); rot.imag = -sinf(coh->phi_[r][c]);
            coh->rx_symb_buf[data_r][c] = cmult(coh->rx_symb_buf[data_r][c], rot);

            /* demodulate */

            i = c*NSYMROW + r;
            rot = cmult(rx_symb[data_r][c], pi_on_4);
            rx_bits[2*i]   = rot.real < 0;
            rx_bits[2*i+1] = rot.imag < 0;

            printf("  c: %d ch_est: %f %f phi_: %f amp_: %f\n",c,  ch_est.real, ch_est.imag, coh->phi_[r][c], coh->amp_[r][c]);
        }
        //exit(0);
    }

    /* shift buffers */

    memcpy(&coh->rx_pilot_buf[0][0], &coh->rx_pilot_buf[NPILOTSFRAME][0], sizeof(COMP)*2*NPILOTSFRAME*PILOTS_NC);
    memcpy(&coh->rx_symb_buf[0][0], &coh->rx_symb_buf[NSYMROW][0], sizeof(COMP)*2*NSYMROW*PILOTS_NC);

}
