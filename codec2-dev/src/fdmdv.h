/*---------------------------------------------------------------------------*\
                                                                             
  FILE........: fdmdv.h
  AUTHOR......: David Rowe                                                          
  DATE CREATED: April 14 2012
                                                                             
  Header file for a Frequency Divison Multiplexed Modem for Digital
  Voice (FDMDV) over HF channels.
                                                                             
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

struct FDMDV;
    
struct FDMDV *fdmdv_create(void);
void          fdmdv_destroy(struct FDMDV *fdmdv_state);
    
void          fdmdv_mod(struct FDMDV *fdmdv_state, int tx_bits[], COMP tx_fdm[]);
void          fdmdv_demod(struct FDMDV *fdmdv_state, int tx_bits[], float rx_fdm[], int *nin);
    
void          fdmdv_get_test_bits(struct FDMDV *fdmdv_state, int tx_bits[]);
void          fdmdv_put_test_bits(struct FDMDV *fdmdv_state, int rx_bits[]);
    
float         fdmdv_get_snr(struct FDMDV *fdmdv_state);
void          fdmdv_get_waterfall_line(struct FDMDV *fdmdv_state, float magnitudes[], int *magnitude_points);

#endif

#ifdef __cplusplus
}
#endif
