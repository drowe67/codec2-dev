/*---------------------------------------------------------------------------*\
                                                                             
  FILE........: fdmdv.h
  AUTHOR......: David Rowe                                                          
  DATE CREATED: April 14 2012
                                                                             
  A 1400 bit/s Frequency Divison Multiplexed Digital Voice (FDMDV)
  modem.  Used for digital audio over HF SSB. The input and output
  sample rate is 8000 Hz. See tfdmdv_mod.c and tfdmdv_demod.c for
  example usage.
                     
  References:
 
    [1] http://n1su.com/fdmdv/FDMDV_Docs_Rel_1_4b.pdf

\*---------------------------------------------------------------------------*/

/*
  Copyright (C) 2012 David Rowe

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

#ifndef __FDMDV__
#define __FDMDV__

#ifdef __cplusplus
extern "C" {
#endif

#include "comp.h"

#define FDMDV_BITS_PER_FRAME      28  /* odd/even frames 56 bits, 1400 bit/s  */
#define FDMDV_SAMPLES_PER_FRAME  160  /* 8000 Hz sample rate                  */
#define FDMDV_SCALE             1000  /* suggested scaling for 16 bit shorts  */
#define FDMDV_NSYM                15

struct FDMDV;
    
struct FDMDV_STATS {
    float  snr;                    /* estimated SNR of rx signal in dB                   */
    COMP   rx_symbols[FDMDV_NSYM]; /* latest received symbols, for scatter plot          */ 
    int    fest_coarse_fine;       /* freq est state, 0-coarse 1-fine                    */ 
    float  foff;                   /* estimated freq offset in Hz                        */       
    float  rx_timing;              /* timing offset -1..1 as fraction of symbol period   */
    float  clock_offset;           /* Estimated tx/rx sample clock offset in ppm         */
};

struct FDMDV *fdmdv_create(void);
void          fdmdv_destroy(struct FDMDV *fdmdv_state);
    
void          fdmdv_mod(struct FDMDV *fdmdv_state, COMP tx_fdm[], int tx_bits[], int *sync_bit);
void          fdmdv_demod(struct FDMDV *fdmdv_state, int rx_bits[], int *sync_bit, float rx_fdm[], int *nin);
    
void          fdmdv_get_test_bits(struct FDMDV *fdmdv_state, int tx_bits[]);
void          fdmdv_put_test_bits(struct FDMDV *f, int *sync, int *bit_errors, int rx_bits[]);
    
void          fdmdv_get_demod_stats(struct FDMDV *fdmdv_state, struct FDMDV_STATS *fdmdv_stats);
void          fdmdv_get_waterfall_line(struct FDMDV *fdmdv_state, float magnitudes[], int *magnitude_points);

#endif

#ifdef __cplusplus
}
#endif
