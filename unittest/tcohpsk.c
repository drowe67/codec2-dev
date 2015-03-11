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
#include "test_bits_coh.h"
#include "octave.h"

#define FRAMES 35
#define CHANNEL_BUF_SIZE (10*M)

int main(int argc, char *argv[])
{
    int           tx_bits[COHPSK_BITS_PER_FRAME];
    COMP          tx_symbols[NSYMROWPILOT][PILOTS_NC];
    
    int           tx_bits_log[COHPSK_BITS_PER_FRAME*FRAMES];
    COMP          tx_symbols_log[NSYMROWPILOT*FRAMES][PILOTS_NC];
                                          
    FILE         *fout;
    int           f, r,c,rx_sym_log_r;

    rx_sym_log_r=0;

    for(f=0; f<FRAMES; f++) {
        
	/* --------------------------------------------------------*\
	                          Modulator
	\*---------------------------------------------------------*/

	bits_to_qpsk_symbols(tx_symbols, (int*)test_bits_coh, sizeof(test_bits_coh));

	/* --------------------------------------------------------*\
	                       Log each vector 
	\*---------------------------------------------------------*/

	memcpy(&tx_bits_log[COHPSK_BITS_PER_FRAME*f], tx_bits, sizeof(int)*COHPSK_BITS_PER_FRAME);
	for(r=0; r<NSYMROWPILOT; r++, rx_sym_log_r++) {
            for(c=0; c<PILOTS_NC; c++) 
		tx_symbols_log[rx_sym_log_r][c] = tx_symbols[r][c]; 
        }
	assert(rx_sym_log_r <= NSYMROWPILOT*FRAMES);
    }


    /*---------------------------------------------------------*\
               Dump logs to Octave file for evaluation 
                      by tcohpsk.m Octave script
    \*---------------------------------------------------------*/

    fout = fopen("tcohpsk_out.txt","wt");
    assert(fout != NULL);
    fprintf(fout, "# Created by tcohpsk.c\n");
    octave_save_int(fout, "tx_bits_log_c", tx_bits_log, 1, COHPSK_BITS_PER_FRAME*FRAMES);
    octave_save_complex(fout, "tx_symbols_log_c", (COMP*)tx_symbols_log, NSYMROWPILOT*FRAMES, PILOTS_NC, PILOTS_NC);  
    fclose(fout);

    return 0;
}

