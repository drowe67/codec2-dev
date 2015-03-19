/*---------------------------------------------------------------------------*\
                                                                             
  FILE........: tcohpsk.c
  AUTHOR......: David Rowe  
  DATE CREATED: March 2015
                                                                             
  Tests for the C version of the coherent PSK FDM modem.  This program
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
#include "comp_prim.h"
//#include "noise_samples.h"

#define FRAMES 35
#define RS     50
#define FOFF   1

int main(int argc, char *argv[])
{
    struct COHPSK *coh;
    int            tx_bits[COHPSK_BITS_PER_FRAME];
    COMP           tx_symb[NSYMROWPILOT][PILOTS_NC];
    COMP           tx_fdm[M*NSYMROWPILOT];
    COMP           rx_fdm[M*NSYMROWPILOT];
    COMP           ch_symb[NSYMROWPILOT][PILOTS_NC];
    int            rx_bits[COHPSK_BITS_PER_FRAME];
    
    int            tx_bits_log[COHPSK_BITS_PER_FRAME*FRAMES];
    COMP           tx_symb_log[NSYMROWPILOT*FRAMES][PILOTS_NC];
    COMP           tx_fdm_log[M*NSYMROWPILOT*FRAMES];
    COMP           rx_fdm_log[M*NSYMROWPILOT*FRAMES];
    COMP           ch_symb_log[NSYMROWPILOT*FRAMES][PILOTS_NC];

    float          rx_amp_log[NSYMROW*FRAMES][PILOTS_NC];
    float          rx_phi_log[NSYMROW*FRAMES][PILOTS_NC];
    COMP           rx_symb_log[NSYMROW*FRAMES][PILOTS_NC];
    int            rx_bits_log[COHPSK_BITS_PER_FRAME*FRAMES];
                                          
    FILE          *fout;
    int            f, r, c, log_r, log_data_r, noise_r, i;
    COMP           phase, freq;
    int           *ptest_bits_coh, *ptest_bits_coh_end;
    float           freq_hz;

    struct FDMDV  *fdmdv;
    COMP           rx_baseband[PILOTS_NC][M+M/P];
    int            nin;
    COMP           rx_filt[PILOTS_NC][P+1];
    COMP           rx_filt_log[PILOTS_NC][(P+1)*NSYMROWPILOT*FRAMES];
    int            rx_filt_log_col_index = 0;
    float          env[NT*P];
    COMP           rx_filter_memory[PILOTS_NC][NFILTER];
    float          rx_timing;
    COMP           tx_onesym[PILOTS_NC];
    COMP           rx_onesym[PILOTS_NC];
    int            rx_baseband_log_col_index = 0;
    COMP           rx_baseband_log[PILOTS_NC][(M+M/P)*NSYMROWPILOT*FRAMES];

    coh = cohpsk_create();
    assert(coh != NULL);

    log_r = log_data_r = noise_r = 0;
    ptest_bits_coh = (int*)test_bits_coh;
    ptest_bits_coh_end = (int*)test_bits_coh + sizeof(test_bits_coh)/sizeof(int);

    memcpy(tx_bits, test_bits_coh, sizeof(int)*COHPSK_BITS_PER_FRAME);

    phase.real = 1.0; phase.imag = 0.0; 
    freq.real = cos(2.0*M_PI*FOFF/RS); freq.imag = sin(2.0*M_PI*FOFF/RS);

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
            rx_filter_memory[c][i].real = 0.0;
            rx_filter_memory[c][i].imag = 0.0;
        }
    }

    /* Main Loop ---------------------------------------------------------------------*/

    for(f=0; f<FRAMES; f++) {
        
	/* --------------------------------------------------------*\
	                          Mod
	\*---------------------------------------------------------*/

        memcpy(tx_bits, ptest_bits_coh, sizeof(int)*COHPSK_BITS_PER_FRAME);
        ptest_bits_coh += COHPSK_BITS_PER_FRAME;
        if (ptest_bits_coh >= ptest_bits_coh_end)
            ptest_bits_coh = (int*)test_bits_coh;
	bits_to_qpsk_symbols(tx_symb, (int*)tx_bits, COHPSK_BITS_PER_FRAME);

        for(r=0; r<NSYMROWPILOT; r++) {
           for(c=0; c<PILOTS_NC; c++) 
              tx_onesym[c] = tx_symb[r][c];         
          tx_filter_and_upconvert(&tx_fdm[r*M], fdmdv->Nc , tx_onesym, fdmdv->tx_filter_memory, 
                                  fdmdv->phase_tx, fdmdv->freq, &fdmdv->fbb_phase_tx, fdmdv->fbb_rect);
        }

        #ifdef TMP
        for(r=0; r<NSYMROWPILOT; r++,noise_r++) {
            phase = cmult(phase,freq);
            for(c=0; c<PILOTS_NC; c++) {
                ch_symb[r][c] = cadd(cmult(tx_symb[r][c], phase), noise[noise_r][c]);
                //printf("%d %d %f %f \n",r,c,ch_symb[r][c].real, ch_symb[r][c].imag);
                //ch_symb[r][c] = cmult(tx_symb[r][c], phase);
             }
        }
        phase = fcmult(1.0/cabsolute(phase), phase);
        #endif

        memcpy(rx_fdm, tx_fdm, sizeof(COMP)*M*NSYMROWPILOT);

	/* --------------------------------------------------------*\
	                          Demod
	\*---------------------------------------------------------*/

        nin = M;
        for (r=0; r<NSYMROWPILOT; r++) {
          fdmdv_freq_shift(&rx_fdm[r*M], &rx_fdm[r*M], -FDMDV_FCENTRE, &fdmdv->fbb_phase_rx, nin);
          fdm_downconvert(rx_baseband, fdmdv->Nc, &rx_fdm[r*M], fdmdv->phase_rx, fdmdv->freq, nin);
          rx_filter(rx_filt, fdmdv->Nc, rx_baseband, rx_filter_memory, nin);
	  rx_timing = rx_est_timing(rx_onesym, fdmdv->Nc, rx_filt, fdmdv->rx_filter_mem_timing, env, nin);
          
          for(c=0; c<PILOTS_NC; c++) {
             ch_symb[r][c] = rx_onesym[c];
          }
          
         for(c=0; c<PILOTS_NC; c++) {       
            for(i=0; i<nin; i++) {
              rx_baseband_log[c][rx_baseband_log_col_index + i] = rx_baseband[c][i]; 
            }
	  }
	  rx_baseband_log_col_index += nin;        

         //if (f == 3)
          //    exit(0);
 	  for(c=0; c<PILOTS_NC; c++) {       
            for(i=0; i<P; i++) {
              rx_filt_log[c][rx_filt_log_col_index + i] = rx_filt[c][i]; 
            }
	  }
	  rx_filt_log_col_index += P;        

        }
        qpsk_symbols_to_bits(coh, rx_bits, ch_symb);

	/* --------------------------------------------------------*\
	                       Log each vector 
	\*---------------------------------------------------------*/

	memcpy(&tx_bits_log[COHPSK_BITS_PER_FRAME*f], tx_bits, sizeof(int)*COHPSK_BITS_PER_FRAME);
	memcpy(&tx_fdm_log[M*NSYMROWPILOT*f], tx_fdm, sizeof(COMP)*M*NSYMROWPILOT);
       	memcpy(&rx_fdm_log[M*NSYMROWPILOT*f], rx_fdm, sizeof(COMP)*M*NSYMROWPILOT);

	for(r=0; r<NSYMROWPILOT; r++, log_r++) {
            for(c=0; c<PILOTS_NC; c++) {
		tx_symb_log[log_r][c] = tx_symb[r][c]; 
		ch_symb_log[log_r][c] = ch_symb[r][c]; 
            }
        }

	for(r=0; r<NSYMROW; r++, log_data_r++) {
            for(c=0; c<PILOTS_NC; c++) {
		rx_amp_log[log_data_r][c] = coh->amp_[r][c]; 
		rx_phi_log[log_data_r][c] = coh->phi_[r][c]; 
		rx_symb_log[log_data_r][c] = coh->rx_symb_buf[r][c]; 
            }
        }
	memcpy(&rx_bits_log[COHPSK_BITS_PER_FRAME*f], rx_bits, sizeof(int)*COHPSK_BITS_PER_FRAME);

	assert(log_r <= NSYMROWPILOT*FRAMES);
	assert(noise_r <= NSYMROWPILOT*FRAMES);
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
    octave_save_complex(fout, "tx_fdm_log_c", (COMP*)tx_fdm_log, 1, M*NSYMROWPILOT*FRAMES, M*NSYMROWPILOT*FRAMES);  
    octave_save_complex(fout, "rx_fdm_log_c", (COMP*)rx_fdm_log, 1, M*NSYMROWPILOT*FRAMES, M*NSYMROWPILOT*FRAMES);  
    octave_save_complex(fout, "rx_baseband_log_c", (COMP*)rx_baseband_log, PILOTS_NC, rx_baseband_log_col_index, (M+M/P)*FRAMES*NSYMROWPILOT);  
    octave_save_complex(fout, "rx_filt_log_c", (COMP*)rx_filt_log, PILOTS_NC, rx_filt_log_col_index, (P+1)*FRAMES*NSYMROWPILOT);  
    octave_save_complex(fout, "ch_symb_log_c", (COMP*)ch_symb_log, NSYMROWPILOT*FRAMES, PILOTS_NC, PILOTS_NC);  
    octave_save_float(fout, "rx_amp_log_c", (float*)rx_amp_log, NSYMROW*FRAMES, PILOTS_NC, PILOTS_NC);  
    octave_save_float(fout, "rx_phi_log_c", (float*)rx_phi_log, NSYMROW*FRAMES, PILOTS_NC, PILOTS_NC);  
    octave_save_complex(fout, "rx_symb_log_c", (COMP*)rx_symb_log, NSYMROW*FRAMES, PILOTS_NC, PILOTS_NC);  
    octave_save_int(fout, "rx_bits_log_c", rx_bits_log, 1, COHPSK_BITS_PER_FRAME*FRAMES);
    fclose(fout);

    fdmdv_destroy(fdmdv);
    cohpsk_destroy(coh);

    return 0;
}

