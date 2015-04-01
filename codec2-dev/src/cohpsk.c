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
#include "fdmdv_internal.h"
#include "pilots_coh.h"
#include "comp_prim.h"
#include "kiss_fft.h"

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

    /* set up buffer of tx pilot symbols for coh demod on rx */

    for(r=0; r<2*NPILOTSFRAME; ) {
        for(p=0; p<NPILOTSFRAME; r++, p++) {
            for(c=0; c<PILOTS_NC; c++) {
                coh->pilot2[r][c] = pilots_coh[p][c];
            }
        }
    }
    
    /* coarse freq offset FFT init */

    coh->fft_coarse_fest = kiss_fft_alloc (COARSE_FEST_NDFT, 0, NULL, NULL);
    assert(coh->fft_coarse_fest != NULL);

    /* Clear symbol buffer memory */

    for (r=0; r<NCT_SYMB_BUF; r++) {
        for(c=0; c<PILOTS_NC; c++) {
            coh->ct_symb_buf[r][c].real = 0.0;
            coh->ct_symb_buf[r][c].imag = 0.0;
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
    KISS_FFT_FREE(coh->fft_coarse_fest);    
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


/*---------------------------------------------------------------------------*\
                                                       
  FUNCTION....: coarse_freq_offset_est()	     
  AUTHOR......: David Rowe			      
  DATE CREATED: March 2015

  Determines an estimate of frequency offset, advances to next sync state.

  TODO: This is currently very stack heavy for an embedded uC.  Tweak
        algorithm in Octave and C to use a smaller DFT and test.

\*---------------------------------------------------------------------------*/

void coarse_freq_offset_est(struct COHPSK *coh, struct FDMDV *fdmdv, COMP ch_fdm_frame[], int sync, int *next_sync)
{
    float f_start, f_stop, sc, h, magsq, num, den, bin_est;
    COMP  s[COARSE_FEST_NDFT], S[COARSE_FEST_NDFT];
    int   bin_start, bin_stop, i;

    if (sync == 0) {
        f_start = FDMDV_FCENTRE - (((float)(PILOTS_NC-1)/2)+2)*FSEP; 
        f_stop = FDMDV_FCENTRE + (((float)(PILOTS_NC-1)/2)+2)*FSEP;
        sc = (float)COARSE_FEST_NDFT/FS;
        bin_start = floorf(f_start*sc+0.5);
        bin_stop = floorf(f_stop*sc+0.5);
        // printf("f_start: %f f_stop: %f sc: %f bin_start: %d bin_stop: %d\n",
        //       f_start, f_stop, sc, bin_start, bin_stop);

        for(i=0; i<NSYMROWPILOT*M; i++) {
            h = 0.5 - 0.5*cosf(2*M_PI*i/(NSYMROWPILOT*M-1));
            s[i] = fcmult(h, ch_fdm_frame[i]); 
        }
        for (; i<COARSE_FEST_NDFT; i++) {
            s[i].real = 0.0;
            s[i].imag = 0.0;
        }

        kiss_fft(coh->fft_coarse_fest, (kiss_fft_cpx *)s, (kiss_fft_cpx *)S);
        
        /* find centroid of signal energy inside search window */

        num = den = 0.0;
        for(i=bin_start; i<bin_stop; i++) {
            magsq = S[i].real*S[i].real + S[i].imag*S[i].imag;
            num += (float)(i+1)*magsq;
            den += magsq;
        }
        bin_est = num/den;
        coh->f_est = bin_est/sc;

        printf("bin_est: %f coarse freq est: %f\n", bin_est, coh->f_est);
     
        *next_sync = 1;
    }
}


/*---------------------------------------------------------------------------*\
                                                       
  FUNCTION....: frame_sync_fine_timing_est()	     
  AUTHOR......: David Rowe			      
  DATE CREATED: April 2015

  Returns an estimate of frequency offset, advances to next sync state.

  TODO: This is very stack heavy for an embedded uC.  Tweak algorthim to use
        a smaller DFT and test.

\*---------------------------------------------------------------------------*/

void frame_sync_fine_timing_est(struct COHPSK *coh, COMP ch_symb[][PILOTS_NC], int sync, int *next_sync)
{
    int   sampling_points[] = {0, 1, 6, 7};
    int   r,c,i,p,t;
    float f_fine, mag, max_corr, max_mag;
    COMP  f_fine_rect, f_corr, corr;

    /* update memory in symbol buffer */

    for (r=0; r<NCT_SYMB_BUF-NSYMROWPILOT; r++) {
        for(c=0; c<PILOTS_NC; c++) {
            coh->ct_symb_buf[r][c] = coh->ct_symb_buf[r+NSYMROWPILOT][c];
        }
    }

    for (r=NCT_SYMB_BUF-NSYMROWPILOT,i=0; r<NCT_SYMB_BUF; r++,i++) {
       for(c=0; c<PILOTS_NC; c++) {
           coh->ct_symb_buf[r][c] = ch_symb[i][c];
           //printf("%d %d %f %f\n", i,c,ch_symb[i][c].real, ch_symb[i][c].imag);
       }
    }

    /* sample pilots at start of this frame and start of next frame */

    if (sync == 2) {

        /* sample correlation over 2D grid of time and fine freq points */

        max_corr = 0;
        for (f_fine=-20; f_fine<=20; f_fine+=1.0) {
            for (t=0; t<NSYMROWPILOT; t++) {
                corr.real = 0.0; corr.imag = 0.0; mag = 0.0;
                for (c=0; c<PILOTS_NC; c++) {
                    for (p=0; p<NPILOTSFRAME+2; p++) {
                        f_fine_rect.real = cosf(-f_fine*2.0*M_PI*(sampling_points[p]+1.0)/RS);
                        f_fine_rect.imag = sinf(-f_fine*2.0*M_PI*(sampling_points[p]+1.0)/RS);
                        f_corr = cmult(f_fine_rect, coh->ct_symb_buf[t+sampling_points[p]][c]);
                        corr = cadd(corr, fcmult(coh->pilot2[p][c], f_corr));
                        mag  += cabsolute(f_corr);
                    }
                }
                //printf("  f: %f  t: %d corr: %f %f\n", f_fine, t, real(corr), imag(corr));
                if (cabsolute(corr) > max_corr) {
                    max_corr = cabsolute(corr);
                    max_mag = mag;
                    coh->ct = t;
                    coh->f_fine_est = f_fine;
                }
            }
        }


        coh->ff_rect.real = cosf(-coh->f_fine_est*2.0*M_PI/RS);
        coh->ff_rect.imag = sinf(-coh->f_fine_est*2.0*M_PI/RS);
        printf("  fine freq f: %f max_corr: %f max_mag: %f ct: %d\n", coh->f_fine_est, max_corr, max_mag, coh->ct);
 
        if (max_corr/max_mag > 0.9) {
            printf("in sync!\n");
            *next_sync = 4;
        }
        else {
            *next_sync = 0;
            printf("  back to coarse freq offset ets...\n");
        }
        //exit(0);
    }
}


int sync_state_machine(sync, next_sync)
{
    if (sync == 1)
        next_sync = 2;
    sync = next_sync;
}



