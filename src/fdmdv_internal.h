/*---------------------------------------------------------------------------*\
                                                                             
  FILE........: fdmdv_internal.h
  AUTHOR......: David Rowe                                                          
  DATE CREATED: April 16 2012
                                                                             
  Header file for FDMDV internal functions, exposed via this header
  file for testing.
                                                                             
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

#ifndef __FDMDV_INTERNAL__
#define __FDMDV_INTERNAL__

#include "comp.h"

/*---------------------------------------------------------------------------*\
                                                                             
                               DEFINES

\*---------------------------------------------------------------------------*/

#define FS                    8000  /* sample rate in Hz                                                    */
#define T                 (1.0/FS)  /* sample period in seconds                                             */
#define RS                      50  /* symbol rate in Hz                                                    */
#define NC                      14  /* number of carriers                                                   */
#define NB                       2  /* Bits/symbol for QPSK modulation                                      */
#define RB              (NC*RS*NB)  /* bit rate                                                             */
#define M                  (FS/RS)  /* oversampling factor                                                  */
#define NSYM                     6  /* number of symbols to filter over                                     */
#define FSEP                    75  /* Separation between carriers (Hz)                                     */
#define FCENTRE               1200  /* Centre frequency, Nc/2 carriers below this, Nc/2 carriers above (Hz) */
#define NT                       5  /* number of symbols we estimate timing over                            */
#define P                        4  /* oversample factor used for initial rx symbol filtering               */
#define NFILTER            (NSYM*M) /* size of tx/rx filters at sample rate M                               */
#define NFILTERTIMING (M+Nfilter+M) /* filter memory used for resampling after timing estimation            */

#define NTEST_BITS        (NC*NB*4) /* length of test bit sequence */

#define PI              3.141592654

/*---------------------------------------------------------------------------*\
                                                                             
                               STRUCT for States

\*---------------------------------------------------------------------------*/

struct FDMDV {
    int  current_test_bit;
    int  tx_pilot_bit;
    COMP prev_tx_symbols[NC+1];
    COMP tx_filter_memory[NC+1][NFILTER];
    COMP phase_tx[NC+1];
    COMP freq[NC+1];
};

/*---------------------------------------------------------------------------*\
                                                                             
                              FUNCTION PROTOTYPES

\*---------------------------------------------------------------------------*/

void bits_to_dqpsk_symbols(COMP tx_symbols[], COMP prev_tx_symbols[], int tx_bits[], int *pilot_bit);
void tx_filter(COMP tx_baseband[NC+1][M], COMP tx_symbols[], COMP tx_filter_memory[NC+1][NFILTER]);
void fdm_upconvert(COMP tx_fdm[], COMP tx_baseband[NC+1][M], COMP phase_tx[], COMP freq_tx[]);

#endif
