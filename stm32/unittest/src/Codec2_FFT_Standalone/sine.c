/*---------------------------------------------------------------------------*\

  FILE........: sine.c
  AUTHOR......: David Rowe
  DATE CREATED: 19/8/2010

  Sinusoidal analysis and synthesis functions.

\*---------------------------------------------------------------------------*/

/*
  Copyright (C) 1990-2010 David Rowe

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

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "defines.h"
#include "sine.h"
#include "kiss_fft.h"

#define HPF_BETA 0.125

/*---------------------------------------------------------------------------*\

				HEADERS

\*---------------------------------------------------------------------------*/

void hs_pitch_refinement(MODEL *model, COMP Sw[], float pmin, float pmax,
			 float pstep);

/*---------------------------------------------------------------------------*\

				FUNCTIONS

\*---------------------------------------------------------------------------*/

C2CONST c2const_create(int Fs, float framelength_s) {
    C2CONST c2const;

    assert((Fs == 8000) || (Fs = 16000));
    c2const.Fs = Fs;
    c2const.n_samp = round(Fs*framelength_s);
    c2const.max_amp = floor(Fs*P_MIN_S/2);
    c2const.p_min = floor(Fs*P_MIN_S);
    c2const.p_max = floor(Fs*P_MAX_S);
    c2const.m_pitch = floor(Fs*M_PITCH_S);
    c2const.Wo_min = TWO_PI/c2const.p_max;
    c2const.Wo_max = TWO_PI/c2const.p_min;

    if (Fs == 8000) {
        c2const.nw = 279;
    } else {
        c2const.nw = 511;  /* actually a bit shorter in time but lets us maintain constant FFT size */
    }

    c2const.tw = Fs*TW_S;

    /*
    fprintf(stderr, "max_amp: %d m_pitch: %d\n", c2const.n_samp, c2const.m_pitch);
    fprintf(stderr, "p_min: %d p_max: %d\n", c2const.p_min, c2const.p_max);
    fprintf(stderr, "Wo_min: %f Wo_max: %f\n", c2const.Wo_min, c2const.Wo_max);
    fprintf(stderr, "nw: %d tw: %d\n", c2const.nw, c2const.tw);
    */

    return c2const;
}

/*---------------------------------------------------------------------------*\

  FUNCTION....: make_analysis_window
  AUTHOR......: David Rowe
  DATE CREATED: 11/5/94

  Init function that generates the time domain analysis window and it's DFT.

\*---------------------------------------------------------------------------*/

void make_analysis_window(C2CONST *c2const, codec2_fft_cfg fft_fwd_cfg, float w[], COMP W[])
{
  float m;
  COMP  wshift[FFT_ENC];
  COMP  temp;
  int   i,j;
  int   m_pitch = c2const->m_pitch;
  int   nw      = c2const->nw;

  /*
     Generate Hamming window centered on M-sample pitch analysis window

  0            M/2           M-1
  |-------------|-------------|
        |-------|-------|
            nw samples

     All our analysis/synthsis is centred on the M/2 sample.
  */

  m = 0.0;
  for(i=0; i<m_pitch/2-nw/2; i++)
    w[i] = 0.0;
  for(i=m_pitch/2-nw/2,j=0; i<m_pitch/2+nw/2; i++,j++) {
    w[i] = 0.5 - 0.5*cosf(TWO_PI*j/(nw-1));
    m += w[i]*w[i];
  }
  for(i=m_pitch/2+nw/2; i<m_pitch; i++)
    w[i] = 0.0;

  /* Normalise - makes freq domain amplitude estimation straight
     forward */

  m = 1.0/sqrtf(m*FFT_ENC);
  for(i=0; i<m_pitch; i++) {
    w[i] *= m;
  }

  /*
     Generate DFT of analysis window, used for later processing.  Note
     we modulo FFT_ENC shift the time domain window w[], this makes the
     imaginary part of the DFT W[] equal to zero as the shifted w[] is
     even about the n=0 time axis if nw is odd.  Having the imag part
     of the DFT W[] makes computation easier.

     0                      FFT_ENC-1
     |-------------------------|

      ----\               /----
           \             /
            \           /          <- shifted version of window w[n]
             \         /
              \       /
               -------

     |---------|     |---------|
       nw/2              nw/2
  */

  for(i=0; i<FFT_ENC; i++) {
    wshift[i].real = 0.0;
    wshift[i].imag = 0.0;
  }
  for(i=0; i<nw/2; i++)
    wshift[i].real = w[i+m_pitch/2];
  for(i=FFT_ENC-nw/2,j=m_pitch/2-nw/2; i<FFT_ENC; i++,j++)
   wshift[i].real = w[j];

  codec2_fft(fft_fwd_cfg, wshift, W);

  /*
      Re-arrange W[] to be symmetrical about FFT_ENC/2.  Makes later
      analysis convenient.

   Before:


     0                 FFT_ENC-1
     |----------|---------|
     __                   _
       \                 /
        \_______________/

   After:

     0                 FFT_ENC-1
     |----------|---------|
               ___
              /   \
     ________/     \_______

  */


  for(i=0; i<FFT_ENC/2; i++) {
    temp.real = W[i].real;
    temp.imag = W[i].imag;
    W[i].real = W[i+FFT_ENC/2].real;
    W[i].imag = W[i+FFT_ENC/2].imag;
    W[i+FFT_ENC/2].real = temp.real;
    W[i+FFT_ENC/2].imag = temp.imag;
  }

}

