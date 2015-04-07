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
    struct FDMDV  *fdmdv;
    int            r,c,p,i;
    float          freq_hz;

    assert(COHPSK_SAMPLES_PER_FRAME == M*NSYMROWPILOT);

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

    coh->ff_phase.real = 1.0; coh->ff_phase.imag = 0.0;
    coh->sync = 0;

    /* set up fdmdv states so we can use those modem functions */

    fdmdv = fdmdv_create(PILOTS_NC - 1);
    for(c=0; c<PILOTS_NC; c++) {
	fdmdv->phase_tx[c].real = 1.0;
 	fdmdv->phase_tx[c].imag = 0.0;

        freq_hz = fdmdv->fsep*( -PILOTS_NC/2 - 0.5 + c + 1.0 );
	fdmdv->freq[c].real = cosf(2.0*M_PI*freq_hz/FS);
 	fdmdv->freq[c].imag = sinf(2.0*M_PI*freq_hz/FS);
 	fdmdv->freq_pol[c]  = 2.0*M_PI*freq_hz/FS;

        //printf("c: %d %f %f\n",c,freq_hz,fdmdv->freq_pol[c]);
        for(i=0; i<NFILTER; i++) {
            coh->rx_filter_memory[c][i].real = 0.0;
            coh->rx_filter_memory[c][i].imag = 0.0;
        }
    }
    coh->fdmdv = fdmdv;

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
    fdmdv_destroy(coh->fdmdv);
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
    int   p, r, c, i;
    COMP  corr, rot, pi_on_4, phi_rect;
    float mag, phi_, amp_;
    short sampling_points[] = {0, 1, 6, 7};

    pi_on_4.real = cosf(M_PI/4); pi_on_4.imag = sinf(M_PI/4);
   
    /* Average pilots to get phase and amplitude estimates we assume
       there are two pilots at the start of each frame and two at the
       end */

    for(c=0; c<PILOTS_NC; c++) {
        corr.real = 0.0; corr.imag = 0.0; mag = 0.0;
        for(p=0; p<NPILOTSFRAME+2; p++) {
            corr = cadd(corr, fcmult(coh->pilot2[p][c], ct_symb_buf[sampling_points[p]][c]));
            mag  += cabsolute(ct_symb_buf[sampling_points[p]][c]);
        }
      
        phi_ = atan2f(corr.imag, corr.real);
        amp_ =  mag/(NPILOTSFRAME+2);
        for(r=0; r<NSYMROW; r++) {
            coh->phi_[r][c] = phi_;
            coh->amp_[r][c] = amp_;
        }
    }

    /* now correct phase of data symbols and make decn on bits */

    for(c=0; c<PILOTS_NC; c++) {
        phi_rect.real = cosf(coh->phi_[0][c]); phi_rect.imag = -sinf(coh->phi_[0][c]);
        //rot.real = 1.0; rot.imag = 0.0;
        for (r=0; r<NSYMROW; r++) {
            i = c*NSYMROW + r;
            coh->rx_symb[r][c] = cmult(ct_symb_buf[NPILOTSFRAME + r][c], phi_rect);
            //printf("%d %d %f %f\n", r,c, coh->rx_symb[r][c].real, coh->rx_symb[r][c].imag);
            //printf("phi_ %d %d %f %f\n", r,c, ct_symb_buf[NPILOTSFRAME + r][c].real, ct_symb_buf[NPILOTSFRAME + r][c].imag);
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
        //printf("f_start: %f f_stop: %f sc: %f bin_start: %d bin_stop: %d\n",
        //        f_start, f_stop, sc, bin_start, bin_stop);

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
        for(i=bin_start; i<=bin_stop; i++) {
            magsq = S[i].real*S[i].real + S[i].imag*S[i].imag;
            num += (float)(i)*magsq;
            den += magsq;
        }
        bin_est = num/den;
        coh->f_est = floor(bin_est/sc+0.5);

        fprintf(stderr, "coarse freq est: %f\n", coh->f_est);
        
        *next_sync = 1;
    }
}


/*---------------------------------------------------------------------------*\
                                                       
  FUNCTION....: frame_sync_fine_timing_est()	     
  AUTHOR......: David Rowe			      
  DATE CREATED: April 2015

  Returns an estimate of frame sync (coarse timing) offset and fine
  frequency offset, advances to next sync state if we have a reliable
  match for frame sync.

  TODO: This is very stack heavy for an embedded uC.  Tweak algorthim to use
        a smaller DFT and test.

\*---------------------------------------------------------------------------*/

void frame_sync_fine_freq_est(struct COHPSK *coh, COMP ch_symb[][PILOTS_NC], int sync, int *next_sync)
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
                        f_fine_rect.real = cosf(f_fine*2.0*M_PI*(sampling_points[p]+1.0)/RS);
                        f_fine_rect.imag = sinf(f_fine*2.0*M_PI*(sampling_points[p]+1.0)/RS);
                        f_corr = cmult(f_fine_rect, coh->ct_symb_buf[t+sampling_points[p]][c]);
                        corr = cadd(corr, fcmult(coh->pilot2[p][c], f_corr));
                        mag  += cabsolute(f_corr);
                    }
                }
                //printf("  f: %f  t: %d corr: %f %f\n", f_fine, t, corr.real, corr.imag);
                if (cabsolute(corr) > max_corr) {
                    max_corr = cabsolute(corr);
                    max_mag = mag;
                    coh->ct = t;
                    coh->f_fine_est = f_fine;
                }
            }
        }


        coh->ff_rect.real = cosf(coh->f_fine_est*2.0*M_PI/RS);
        coh->ff_rect.imag = -sinf(coh->f_fine_est*2.0*M_PI/RS);
        fprintf(stderr, "  fine freq f: %f max_corr: %f max_mag: %f ct: %d\n", coh->f_fine_est, max_corr, max_mag, coh->ct);
 
        if (max_corr/max_mag > 0.9) {
            fprintf(stderr, "in sync!\n");
            *next_sync = 4;
        }
        else {
            *next_sync = 0;
            fprintf(stderr, "  back to coarse freq offset est...\n");
        }
        
    }
}


/*---------------------------------------------------------------------------*\
                                                       
  FUNCTION....: fine_freq_correct()	     
  AUTHOR......: David Rowe			      
  DATE CREATED: April 2015

  Fine frequency correction of symbols at rate Rs.

\*---------------------------------------------------------------------------*/

void fine_freq_correct(struct COHPSK *coh, int sync, int next_sync) {
    int   r,c;
    float mag;

  /*
    We can decode first frame that we achieve sync.  Need to fine freq
    correct all of it's symbols, including pilots.  From then on, just
    correct new symbols into frame.  make copy, so if we lose sync we
    havent fine freq corrected ct_symb_buf if next_sync == 4 correct
    all 8 if sync == 2 correct latest 6.
  */

  if ((next_sync == 4) || (sync == 4)) {

      if ((next_sync == 4) && (sync == 2)) {
          
          /* first frame, we've just gotten sync so fine freq correct all Nsymbrowpilot+2 samples */

          for(r=0; r<NSYMROWPILOT+2; r++) {
              coh->ff_phase = cmult(coh->ff_phase, cconj(coh->ff_rect));
              for(c=0; c<PILOTS_NC; c++) {
                  coh->ct_symb_ff_buf[r][c] = coh->ct_symb_buf[coh->ct+r][c];
                  coh->ct_symb_ff_buf[r][c] = cmult(coh->ct_symb_ff_buf[r][c], coh->ff_phase);
              }
          }
      }
      else {
          
          /* second and subsequent frames, just fine freq correct the latest Nsymbrowpilot */

          for(r=0; r<2; r++) {
              for(c=0; c<PILOTS_NC; c++) {
                  coh->ct_symb_ff_buf[r][c] = coh->ct_symb_ff_buf[r+NSYMROWPILOT][c];
              }
          }
          
          for(; r<NSYMROWPILOT+2; r++) {
              coh->ff_phase = cmult(coh->ff_phase, cconj(coh->ff_rect));
              for(c=0; c<PILOTS_NC; c++) {
                  coh->ct_symb_ff_buf[r][c] = coh->ct_symb_buf[coh->ct+r][c];
                  //printf("%d %d %f %f\n", r,c,coh->ct_symb_ff_buf[r][c].real, coh->ct_symb_ff_buf[r][c].imag);
                  coh->ct_symb_ff_buf[r][c] = cmult(coh->ct_symb_ff_buf[r][c], coh->ff_phase);
              }
          }
      }

      mag = cabsolute(coh->ff_phase);
      coh->ff_phase.real /= mag;
      coh->ff_phase.imag /= mag;
  }
}


int sync_state_machine(int sync, int next_sync)
{
    if (sync == 1)
        next_sync = 2;
    sync = next_sync;

    return sync;
}


/*---------------------------------------------------------------------------*\
                                                       
  FUNCTION....: cohpsk_mod()	     
  AUTHOR......: David Rowe			      
  DATE CREATED: 5/4/2015

  COHPSK modulator, take a frame of COHPSK_BITS_PER_FRAME bits and
  generates a frame of COHPSK_SAMPLES_PER_FRAME modulated symbols.

  The output signal is complex to support single sided frequency
  shifting, for example when testing frequency offsets in channel
  simulation.

\*---------------------------------------------------------------------------*/

void cohpsk_mod(struct COHPSK *coh, COMP tx_fdm[], int tx_bits[])
{
    struct FDMDV *fdmdv = coh->fdmdv;
    COMP tx_symb[NSYMROWPILOT][PILOTS_NC];
    COMP tx_onesym[PILOTS_NC];
    int  r,c;

    bits_to_qpsk_symbols(tx_symb, tx_bits, COHPSK_BITS_PER_FRAME);

    for(r=0; r<NSYMROWPILOT; r++) {
        for(c=0; c<PILOTS_NC; c++) 
            tx_onesym[c] = tx_symb[r][c];         
        tx_filter_and_upconvert(&tx_fdm[r*M], fdmdv->Nc , tx_onesym, fdmdv->tx_filter_memory, 
                                fdmdv->phase_tx, fdmdv->freq, &fdmdv->fbb_phase_tx, fdmdv->fbb_rect);
    }
}


/*---------------------------------------------------------------------------*\
                                                       
  FUNCTION....: cohpsk_demod()	     
  AUTHOR......: David Rowe			      
  DATE CREATED: 5/4/2015

  COHPSK demodulator, takes an array of COHPSK_SAMPLES_PER_FRAME
  modulated samples, returns an array of COHPSK_BITS_PER_FRAME bits.

  The input signal is complex to support single sided frequency shifting
  before the demod input (e.g. click to tune feature).

  The number of input samples is fixed, and unlike the FDMDV modem
  doesn't change to adjust for differences in transmit and receive
  sample clocks.  This means frame sync will occasionally be lost,
  however this is hardly noticable for digital voice applications.

  TODO: logic to check if we are still in sync, ride thru bad frames

\*---------------------------------------------------------------------------*/

void cohpsk_demod(struct COHPSK *coh, int rx_bits[], int *reliable_sync_bit, COMP rx_fdm[])
{
    struct FDMDV *fdmdv = coh->fdmdv;
    COMP  rx_fdm_frame_bb[M*NSYMROWPILOT];
    COMP  rx_baseband[PILOTS_NC][M+M/P];
    COMP  rx_filt[PILOTS_NC][P+1];
    float env[NT*P], rx_timing;
    COMP  ch_symb[NSYMROWPILOT][PILOTS_NC];
    COMP  rx_onesym[PILOTS_NC];
    int   sync, next_sync, nin, r, c;

    next_sync = sync = coh->sync;

    coarse_freq_offset_est(coh, fdmdv, rx_fdm, sync, &next_sync);

    /* sample rate demod processing */

    nin = M;
    for (r=0; r<NSYMROWPILOT; r++) {
        fdmdv_freq_shift(&rx_fdm_frame_bb[r*M], &rx_fdm[r*M], -coh->f_est, &fdmdv->fbb_phase_rx, nin);
        fdm_downconvert(rx_baseband, fdmdv->Nc, &rx_fdm_frame_bb[r*M], fdmdv->phase_rx, fdmdv->freq, nin);
        rx_filter(rx_filt, fdmdv->Nc, rx_baseband, coh->rx_filter_memory, nin);
        rx_timing = rx_est_timing(rx_onesym, fdmdv->Nc, rx_filt, fdmdv->rx_filter_mem_timing, env, nin);
          
        for(c=0; c<PILOTS_NC; c++) {
            ch_symb[r][c] = rx_onesym[c];
        }
    }
     
    /* coarse timing (frame sync) and initial fine freq est */
  
    frame_sync_fine_freq_est(coh, ch_symb, sync, &next_sync);
    fine_freq_correct(coh, sync, next_sync);
        
    *reliable_sync_bit = 0;
    if ((sync == 4) || (next_sync == 4)) {
        qpsk_symbols_to_bits(coh, rx_bits, coh->ct_symb_ff_buf);
        *reliable_sync_bit = 1;
    }

    sync = sync_state_machine(sync, next_sync);        

    coh->sync = sync;
}
