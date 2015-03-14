/*---------------------------------------------------------------------------*\
                                                                             
  FILE........: tcopskv.c
  AUTHOR......: David Rowe  
  DATE CREATED: March 2015
                                                                             
  Tests for the C version of the cohernet PSK FDM modem.  This program
  outputs a file of Octave vectors that are loaded and automatically
  tested against the Octave version of the modem by the Octave script
  tcohpsk.m
                                                                             
\*---------------------------------------------------------------------------*/

/*
  Copyright (C) 2015 David Rowe

  All rights reserved.

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU Lesser General Public License version 2, as
  published by the Free Software Foundation.  This program is
  distributed in the hope that it will be useful, but WITHOUT ANY
  WARRANTY; without even the implied warranty of MERCHANTABILITY or
  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public
  License for more details.

  You should have received a copy of the GNU Lesser General Public License
  along with this program; if not, see <http://www.gnu.org/licenses/>.
*/

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "fdmdv_internal.h"
#include "codec2_fdmdv.h"
#include "codec2_cohpsk.h"
#include "cohpsk_defs.h"
#include "cohpsk_internal.h"
#include "test_bits_coh.h"
#include "octave.h"

#define FRAMES 35
#define CHANNEL_BUF_SIZE (10*M)

int main(int argc, char *argv[])
{
    struct COHPSK *coh;
    int            tx_bits[COHPSK_BITS_PER_FRAME];
    COMP           tx_symb[NSYMROWPILOT][PILOTS_NC];
    int            rx_bits[COHPSK_BITS_PER_FRAME];
    
    int            tx_bits_log[COHPSK_BITS_PER_FRAME*FRAMES];
    COMP           tx_symb_log[NSYMROWPILOT*FRAMES][PILOTS_NC];

    float          rx_amp_log[NSYMROW*FRAMES][PILOTS_NC];
    float          rx_phi_log[NSYMROW*FRAMES][PILOTS_NC];
    int            rx_bits_log[COHPSK_BITS_PER_FRAME*FRAMES];
                                          
    FILE          *fout;
    int            f, r, c, log_r, log_data_r;

    coh = cohpsk_create();
    assert(coh != NULL);

    log_r = log_data_r= 0;

    memcpy(tx_bits, test_bits_coh, sizeof(int)*COHPSK_BITS_PER_FRAME);

    for(f=0; f<FRAMES; f++) {
        
	/* --------------------------------------------------------*\
	                          Modem
	\*---------------------------------------------------------*/

	bits_to_qpsk_symbols(tx_symb, (int*)tx_bits, COHPSK_BITS_PER_FRAME);
        qpsk_symbols_to_bits(coh, rx_bits, tx_symb);

	/* --------------------------------------------------------*\
	                       Log each vector 
	\*---------------------------------------------------------*/

	memcpy(&tx_bits_log[COHPSK_BITS_PER_FRAME*f], tx_bits, sizeof(int)*COHPSK_BITS_PER_FRAME);
	for(r=0; r<NSYMROWPILOT; r++, log_r++) {
            for(c=0; c<PILOTS_NC; c++) 
		tx_symb_log[log_r][c] = tx_symb[r][c]; 
        }

	for(r=0; r<NSYMROW; r++, log_data_r++) {
            for(c=0; c<PILOTS_NC; c++) {
		rx_amp_log[log_data_r][c] = coh->amp_[r][c]; 
		rx_phi_log[log_data_r][c] = coh->phi_[r][c]; 
            }
        }
	memcpy(&rx_bits_log[COHPSK_BITS_PER_FRAME*f], rx_bits, sizeof(int)*COHPSK_BITS_PER_FRAME);

	assert(log_r <= NSYMROWPILOT*FRAMES);
	assert(log_data_r <= NSYMROW*FRAMES);
    }

    /*---------------------------------------------------------*\
               Dump logs to Octave file for evaluation 
                      by tcohpsk.m Octave script
    \*---------------------------------------------------------*/

    fout = fopen("tcohpsk_out.txt","wt");
    assert(fout != NULL);
    fprintf(fout, "# Created by tcohpsk.c\n");
    octave_save_int(fout, "tx_bits_log_c", tx_bits_log, 1, COHPSK_BITS_PER_FRAME*FRAMES);
    octave_save_complex(fout, "tx_symb_log_c", (COMP*)tx_symb_log, NSYMROWPILOT*FRAMES, PILOTS_NC, PILOTS_NC);  
    octave_save_complex(fout, "rx_amp_log_c", (COMP*)rx_amp_log, NSYMROW*FRAMES, PILOTS_NC, PILOTS_NC);  
    octave_save_complex(fout, "rx_phi_log_c", (COMP*)rx_phi_log, NSYMROW*FRAMES, PILOTS_NC, PILOTS_NC);  
    fclose(fout);

    cohpsk_destroy(coh);

    return 0;
}

