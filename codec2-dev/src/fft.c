/*---------------------------------------------------------------------------*\
                                                 
  FILE........: fft.c                                                  
  AUTHOR......: Bruce Robertson                                      
  DATE CREATED: 20/11/2010                            
                                                         
  Bridging function to the kiss_fft package.      
                                                               
\*---------------------------------------------------------------------------*/

/*
  Copyright (C) 2010 Bruce Robertson

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

#include <assert.h>
#include "kiss_fft.h"

/*---------------------------------------------------------------------------*\
                                                                            
                                GLOBALS                                       
                                                                             
\*---------------------------------------------------------------------------*/

kiss_fft_cpx *fin        = 0;
kiss_fft_cpx *fout       = 0;
kiss_fft_cfg cfg_forward = 0;
kiss_fft_cfg cfg_reverse = 0;

/*---------------------------------------------------------------------------*\
                                                                             
  initialize_fft(int n)                                                                  
                                                                             
  Initialisation function for kiss_fft. This assumes that all calls to fft() 
  use the same datatypes and are one arrays of the same size.

\*---------------------------------------------------------------------------*/

void
initialize_fft (int n)
{
  assert(!fin && !fout && !cfg_forward && !cfg_reverse);
  fin = KISS_FFT_MALLOC (n * sizeof (kiss_fft_cpx));
  assert(fin != NULL);
  fout = KISS_FFT_MALLOC (n * sizeof (kiss_fft_cpx));
  assert(fout != NULL);
  cfg_forward = kiss_fft_alloc (n, 0, NULL, NULL);
  assert(cfg_forward != NULL);
  cfg_reverse = kiss_fft_alloc (n, 1, NULL, NULL);
  assert(cfg_reverse != NULL);
}

/*---------------------------------------------------------------------------*\
                                                                             
  fft(float x[], int n, int isign)                                                
  Function that calls kiss_fft with the signature of four1 from NRC.

\*---------------------------------------------------------------------------*/


#ifdef BRUCE

/* Efficient but runs into problems if we have two different size FFTs
   in the same program - see notes below */

void
fft (float x[], int n, int isign)
{
  int c;
  kiss_fft_cfg cfg;
  if (cfg_forward == NULL)
    {
      initialize_fft (n);
    }
  for (c = 0; c < n * 2; c += 2)
    {
      fin[c / 2].r = x[c];
      fin[c / 2].i = -x[c + 1];
    }
  if (isign == -1)
    {
      cfg = cfg_reverse;
    }
  else
    {
      cfg = cfg_forward;
    }
  kiss_fft (cfg, fin, fout);
  for (c = 0; c < n * 2; c += 2)
    {
      x[c] = fout[(c) / 2].r;
      x[c + 1] = -fout[(c) / 2].i;
    }
}
#endif

/* This version not as efficient but can handle different size FFTs in
   the same program.  This is reqd in fdmdv and if we link fdmdv and
   codec 2 into same program. If CPU load becomes an issue we could always
   modify to allocate FFT cfg states at start up.

   Or maybe we should just bite the bullet and modify all fft() calls
   to match the kiss_fft calling conventions.  This would mean
   allocating states for each fft at the start of the program which is
   no biggie.

*/

#define DAVID
#ifdef DAVID
void
fft (float x[], int n, int isign)
{
  int             c;
  kiss_fft_cfg    cfg;
  kiss_fft_cpx   *input, *output;

  input = KISS_FFT_MALLOC (n * sizeof (kiss_fft_cpx));
  assert(input != NULL);
  output = KISS_FFT_MALLOC (n * sizeof (kiss_fft_cpx));
  assert(output != NULL);

  for (c = 0; c < n * 2; c += 2) {
      input[c / 2].r = x[c];
      input[c / 2].i = -x[c + 1];
  }
  if (isign == -1)
      cfg = kiss_fft_alloc (n, 1, NULL, NULL);
  else
      cfg = kiss_fft_alloc (n, 0, NULL, NULL);
  kiss_fft (cfg, input, output);
  for (c = 0; c < n * 2; c += 2) {
      x[c] = output[(c) / 2].r;
      x[c + 1] = -output[(c) / 2].i;
  }
  KISS_FFT_FREE(input);
  KISS_FFT_FREE(output);
  KISS_FFT_FREE(cfg);
}
#endif

void cleanup_fft(void)
{
    KISS_FFT_FREE(fin);
    KISS_FFT_FREE(fout);
    KISS_FFT_FREE(cfg_forward);
    KISS_FFT_FREE(cfg_reverse);
}
