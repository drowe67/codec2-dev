/*---------------------------------------------------------------------------*\
                                                                             
  FILE........: fm.c
  AUTHOR......: David Rowe
  DATE CREATED: February 2015
                                                                             
  Functions that implement analog FM modulation and demodulation, see
  also octave/fm.m.
                                                                             
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
                                                                             
                               DEFINES

\*---------------------------------------------------------------------------*/

#define FILT_MEM 200

/*---------------------------------------------------------------------------*\
                                                                             
                               INCLUDES

\*---------------------------------------------------------------------------*/

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "codec2_fm.h"
#include "fm_fir_coeff.h"

/*---------------------------------------------------------------------------*\
                                                                             
                               FUNCTIONS

\*---------------------------------------------------------------------------*/

static COMP cconj(COMP a)
{
    COMP res;

    res.real = a.real;
    res.imag = -a.imag;

    return res;
}

static COMP cmult(COMP a, COMP b)
{
    COMP res;

    res.real = a.real*b.real - a.imag*b.imag;
    res.imag = a.real*b.imag + a.imag*b.real;

    return res;
}

static COMP fcmult(float a, COMP b)
{
    COMP res;

    res.real = a*b.real;
    res.imag = a*b.imag;

    return res;
}

static float cabsolute(COMP a)
{
    return sqrtf(powf(a.real, 2.0) + powf(a.imag, 2.0));
}

/*---------------------------------------------------------------------------*\
                                                       
  FUNCTION....: fm_create	     
  AUTHOR......: David Rowe			      
  DATE CREATED: 24 Feb 2015

  Create and initialise an instance of the "modem".  Returns a pointer
  to the modem states or NULL on failure.  One set of states is
  sufficient for a full duplex modem.

\*---------------------------------------------------------------------------*/

struct FM *fm_create(int nsam)
{
    struct FM *fm;

    fm = (struct FM*)malloc(sizeof(struct FM));
    if (fm == NULL)
	return NULL;
    fm->rx_bb = (COMP*)malloc(sizeof(COMP)*(FILT_MEM+nsam));
    assert(fm->rx_bb != NULL);

    fm->rx_bb_filt_prev.real = 0.0;
    fm->rx_bb_filt_prev.imag = 0.0;
    fm->lo_phase.real = 1.0;
    fm->lo_phase.imag = 0.0;

    fm->rx_dem_mem = (float*)malloc(sizeof(float)*(FILT_MEM+nsam));
    assert(fm->rx_dem_mem != NULL);

    fm->nsam = nsam;

    return fm;
}


void fm_destroy(struct FM *fm_states)
{
    free(fm_states->rx_bb);
    free(fm_states->rx_dem_mem);
    free(fm_states);
}

/*---------------------------------------------------------------------------*\
                                                       
  FUNCTION....: fm_demod	     
  AUTHOR......: David Rowe			      
  DATE CREATED: 24 Feb 2015

  Demodulate a FM signal to baseband audio.

\*---------------------------------------------------------------------------*/

void fm_demod(struct FM *fm_states, float rx_out[], float rx[])
{
  float  Fs = fm_states->Fs;
  float  fc = fm_states->fc; 
  float  wc = 2*M_PI*fc/Fs;
  float  fd = fm_states->fd; 
  float  wd = 2*M_PI*fd/Fs;
  COMP  *rx_bb = fm_states->rx_bb + FILT_MEM;
  COMP   wc_rect, rx_bb_filt, rx_bb_diff;
  float  rx_dem, acc;
  float *rx_dem_mem = fm_states->rx_dem_mem + FILT_MEM;
  int    nsam = fm_states->nsam;
  float  mag;
  int    i,k;

  wc_rect.real = cos(wc); wc_rect.imag = -sin(wc);

  for(i=0; i<nsam; i++) {

      /* down to complex baseband */

      fm_states->lo_phase = cmult(fm_states->lo_phase, wc_rect);
      rx_bb[i] = fcmult(rx[i], fm_states->lo_phase);      

      /* input FIR filter */

      rx_bb_filt.real = 0.0; rx_bb_filt.imag = 0.0;
      for(k=0; k<FILT_MEM; k++) {
          rx_bb_filt.real += rx_bb[i-k].real * bin[k];
          rx_bb_filt.imag += rx_bb[i-k].imag * bin[k];
      }
      //printf("%f %f %f\n", rx[i], wc_rect.real, wc_rect.imag);
      //printf("%f %f %f\n", rx[i], fm_states->lo_phase.real, fm_states->lo_phase.imag);
      //printf("%f %f %f\n", rx[i], rx_bb[i].real, rx_bb[i].imag);
      //printf("%f %f\n", rx_bb_filt.real, rx_bb_filt.imag);
      /* 
         Differentiate first, in rect domain, then find angle, this
         puts signal on the positive side of the real axis and helps
         atan2() behaive.
      */

      rx_bb_diff = cmult(rx_bb_filt, cconj(fm_states->rx_bb_filt_prev));
      fm_states->rx_bb_filt_prev = rx_bb_filt;
  
      rx_dem = atan2(rx_bb_diff.imag, rx_bb_diff.real);

      /* limit maximum phase jumps, to remove static type noise at low SNRs */

      if (rx_dem > wd)
          rx_dem = wd;
      if (rx_dem < -wd)
          rx_dem = -wd;

      rx_dem *= (1/wd);
      //printf("%f %f\n", rx_bb_diff.real, rx_bb_diff.imag);
      rx_dem_mem[i] = rx_dem;
      acc = 0;
      for(k=0; k<FILT_MEM; k++) {
          acc += rx_dem_mem[i-k] * bout[k];
      }
      rx_out[i] = acc;
  }

  /* update filter memories */

  rx_bb      -= FILT_MEM;
  rx_dem_mem -= FILT_MEM;
  for(i=0; i<FILT_MEM; i++) {
      rx_bb[i] = rx_bb[i+nsam];
      rx_dem_mem[i] = rx_dem_mem[i+nsam];
  }

  /* normalise digital oscillator as the magnitude can drift over time */

  mag = cabsolute(fm_states->lo_phase);
  fm_states->lo_phase.real /= mag;	
  fm_states->lo_phase.imag /= mag;	

}
