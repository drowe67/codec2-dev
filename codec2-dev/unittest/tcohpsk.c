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
#include "noise_samples.h"

#define FRAMES      35                    /* LOG_FRAMES is #defined in cohpsk_internal.h                        */
#define SYNC_FRAMES 12                    /* sync state uses up extra log storage as we reprocess several times */
#define FRAMESL     (SYNC_FRAMES*FRAMES)  /* worst case is every frame is out of sync                           */

#define FOFF        58.7
#define DFOFF       (-0.5/(float)COHPSK_FS)
#define ESNODB      8
#define PPM         -1000

extern float pilots_coh[][PILOTS_NC];

int main(int argc, char *argv[])
{
    struct COHPSK *coh;
    int            tx_bits[COHPSK_BITS_PER_FRAME];
    COMP           tx_symb[NSYMROWPILOT][COHPSK_NC*ND];
    COMP           tx_fdm_frame[COHPSK_M*NSYMROWPILOT];
    COMP           ch_fdm_frame[COHPSK_M*NSYMROWPILOT];
    //COMP           rx_fdm_frame_bb[M*NSYMROWPILOT];
    //COMP           ch_symb[NSYMROWPILOT][COHPSK_NC*ND];
    int            rx_bits[COHPSK_BITS_PER_FRAME];
    
    int            tx_bits_log[COHPSK_BITS_PER_FRAME*FRAMES];
    COMP           tx_symb_log[NSYMROWPILOT*FRAMES][COHPSK_NC*ND];
    COMP           tx_fdm_frame_log[COHPSK_M*NSYMROWPILOT*FRAMES];
    COMP           ch_fdm_frame_log[COHPSK_M*NSYMROWPILOT*FRAMES];
    COMP           ch_fdm_frame_log_out[(COHPSK_M*NSYMROWPILOT+1)*FRAMES];
    //COMP           rx_fdm_frame_bb_log[M*NSYMROWPILOT*FRAMES];
    //COMP           ch_symb_log[NSYMROWPILOT*FRAMES][COHPSK_NC*ND];
    COMP           ct_symb_ff_log[NSYMROWPILOT*FRAMES][COHPSK_NC*ND];
    float          rx_amp_log[NSYMROW*FRAMES][COHPSK_NC*ND];
    float          rx_phi_log[NSYMROW*FRAMES][COHPSK_NC*ND];
    COMP           rx_symb_log[NSYMROW*FRAMES][COHPSK_NC*ND];
    int            rx_bits_log[COHPSK_BITS_PER_FRAME*FRAMES];
                                          
    FILE          *fout;
    int            f, r, c, log_r, log_data_r, noise_r, ff_log_r;
    int           *ptest_bits_coh, *ptest_bits_coh_end;
    double         foff;
    COMP           foff_rect, phase_ch;

    struct FDMDV  *fdmdv;
    //COMP           rx_filt[COHPSK_NC*ND][P+1];
    //int            rx_filt_log_col_index = 0;
    //float          env[NT*P];
    //float           __attribute__((unused)) rx_timing;
    COMP           tx_onesym[COHPSK_NC*ND];
    //COMP           rx_onesym[COHPSK_NC*ND];
    //int            rx_baseband_log_col_index = 0;
    //COMP           rx_baseband_log[COHPSK_NC*ND][(M+M/P)*NSYMROWPILOT*FRAMES];
    float            f_est_log[FRAMES];
    int              f_est_samples;

    int            log_bits;
    float          EsNo, variance;
    COMP           scaled_noise;
    int            reliable_sync_bit;
    int            ch_fdm_frame_log_index, nin_frame, tmp, nout;

    coh = cohpsk_create();
    fdmdv = coh->fdmdv;
    assert(coh != NULL);

    /* these puppies are used for logging data in the bowels on the modem */

    coh->rx_baseband_log_col_sz = (COHPSK_M+COHPSK_M/P)*NSYMROWPILOT*FRAMESL;
    coh->rx_baseband_log = (COMP *)malloc(sizeof(COMP)*COHPSK_NC*ND*coh->rx_baseband_log_col_sz);

    coh->rx_filt_log_col_sz = (P+1)*NSYMROWPILOT*FRAMESL;
    coh->rx_filt_log = (COMP *)malloc(sizeof(COMP)*COHPSK_NC*ND*coh->rx_filt_log_col_sz);

    coh->ch_symb_log_col_sz = COHPSK_NC*ND;
    coh->ch_symb_log = (COMP *)malloc(sizeof(COMP)*NSYMROWPILOT*FRAMESL*coh->ch_symb_log_col_sz);

    coh->rx_timing_log = (float*)malloc(sizeof(float)*NSYMROWPILOT*FRAMESL);
    
    /* init stuff */

    log_r = log_data_r = noise_r = log_bits = ff_log_r = f_est_samples = 0;
    ptest_bits_coh = (int*)test_bits_coh;
    ptest_bits_coh_end = (int*)test_bits_coh + sizeof(test_bits_coh)/sizeof(int);
    memcpy(tx_bits, test_bits_coh, sizeof(int)*COHPSK_BITS_PER_FRAME);

    phase_ch.real = 1.0; phase_ch.imag = 0.0; 
    foff = FOFF;
     
    /*  each carrier has power = 2, total power 2Nc, total symbol rate
        NcRs, noise BW B=Fs Es/No = (C/Rs)/(N/B), N = var =
        2NcFs/NcRs(Es/No) = 2Fs/Rs(Es/No) */

    EsNo = pow(10.0, ESNODB/10.0);
    variance = 2.0*COHPSK_FS/(COHPSK_RS*EsNo);
    //fprintf(stderr, "doff: %e\n", DFOFF);

    /* Main Loop ---------------------------------------------------------------------*/

    for(f=0; f<FRAMES; f++) {

	/* --------------------------------------------------------*\
	                          Mod
	\*---------------------------------------------------------*/

        memcpy(tx_bits, ptest_bits_coh, sizeof(int)*COHPSK_BITS_PER_FRAME);
        ptest_bits_coh += COHPSK_BITS_PER_FRAME;
        if (ptest_bits_coh >= ptest_bits_coh_end) {
            ptest_bits_coh = (int*)test_bits_coh;
        }

	bits_to_qpsk_symbols(tx_symb, (int*)tx_bits, COHPSK_BITS_PER_FRAME);

        for(r=0; r<NSYMROWPILOT; r++) {
           for(c=0; c<COHPSK_NC*ND; c++) 
               tx_onesym[c] = tx_symb[r][c];         
           tx_filter_and_upconvert_coh(&tx_fdm_frame[r*COHPSK_M], COHPSK_NC*ND , tx_onesym, fdmdv->tx_filter_memory, 
                                        fdmdv->phase_tx, fdmdv->freq, &fdmdv->fbb_phase_tx, fdmdv->fbb_rect);
        }
        cohpsk_clip(tx_fdm_frame);

	/* --------------------------------------------------------*\
	                          Channel
	\*---------------------------------------------------------*/


        for(r=0; r<NSYMROWPILOT*COHPSK_M; r++) {
            foff_rect.real = cos(2.0*M_PI*foff/COHPSK_FS); foff_rect.imag = sin(2.0*M_PI*foff/COHPSK_FS);
            foff += DFOFF;
            phase_ch = cmult(phase_ch, foff_rect);
            ch_fdm_frame[r] = cmult(tx_fdm_frame[r], phase_ch);
        }
        phase_ch.real /= cabsolute(phase_ch);
        phase_ch.imag /= cabsolute(phase_ch);
        //fprintf(stderr, "%f\n", foff);
        for(r=0; r<NSYMROWPILOT*COHPSK_M; r++,noise_r++) {
            scaled_noise = fcmult(sqrt(variance), noise[noise_r]);
            ch_fdm_frame[r] = cadd(ch_fdm_frame[r], scaled_noise);
        }

        /* tx vector logging */

	memcpy(&tx_bits_log[COHPSK_BITS_PER_FRAME*f], tx_bits, sizeof(int)*COHPSK_BITS_PER_FRAME);
	memcpy(&tx_fdm_frame_log[COHPSK_M*NSYMROWPILOT*f], tx_fdm_frame, sizeof(COMP)*COHPSK_M*NSYMROWPILOT);
	memcpy(&ch_fdm_frame_log[COHPSK_M*NSYMROWPILOT*f], ch_fdm_frame, sizeof(COMP)*COHPSK_M*NSYMROWPILOT);
  
	for(r=0; r<NSYMROWPILOT; r++, log_r++) {
            for(c=0; c<COHPSK_NC*ND; c++) {
		tx_symb_log[log_r][c] = tx_symb[r][c]; 
            }
        }
    }

    /* Fs offset simulation */

    nout = cohpsk_fs_offset(ch_fdm_frame_log_out, ch_fdm_frame_log, COHPSK_M*NSYMROWPILOT*FRAMES, PPM);
    assert(nout < (COHPSK_M*NSYMROWPILOT+1)*FRAMES);

    nin_frame = COHPSK_SAMPLES_PER_FRAME;
    ch_fdm_frame_log_index = 0;

    for(f=0; f<FRAMES; f++) {
        coh->frame = f;

        //printf("nin_frame: %d\n", nin_frame);

	/* --------------------------------------------------------*\
	                          Demod
	\*---------------------------------------------------------*/

        assert(ch_fdm_frame_log_index < COHPSK_M*NSYMROWPILOT*FRAMES);
        tmp = nin_frame;
        cohpsk_demod(coh, rx_bits, &reliable_sync_bit, &ch_fdm_frame_log_out[ch_fdm_frame_log_index], &nin_frame);
        ch_fdm_frame_log_index += tmp;

 	/* --------------------------------------------------------*\
	                       Log each vector 
	\*---------------------------------------------------------*/

        if (coh->sync == 1) {

            for(r=0; r<NSYMROWPILOT; r++, ff_log_r++) {
                for(c=0; c<COHPSK_NC*ND; c++) {
                    ct_symb_ff_log[ff_log_r][c] = coh->ct_symb_ff_buf[r][c]; 
                }
            }

            for(r=0; r<NSYMROW; r++, log_data_r++) {
                for(c=0; c<COHPSK_NC*ND; c++) {
                    rx_amp_log[log_data_r][c] = coh->amp_[r][c]; 
                    rx_phi_log[log_data_r][c] = coh->phi_[r][c]; 
                    rx_symb_log[log_data_r][c] = coh->rx_symb[r][c]; 
                }
            }
            memcpy(&rx_bits_log[COHPSK_BITS_PER_FRAME*log_bits], rx_bits, sizeof(int)*COHPSK_BITS_PER_FRAME);
            log_bits++;
            f_est_log[f_est_samples++] = coh->f_est;
        }

	assert(log_r <= NSYMROWPILOT*FRAMES);
	assert(noise_r <= NSYMROWPILOT*COHPSK_M*FRAMES);
	assert(log_data_r <= NSYMROW*FRAMES);

        printf("\r  [%d]", f+1);
    }
    printf("\n");

    /*---------------------------------------------------------*\
               Dump logs to Octave file for evaluation 
                      by tcohpsk.m Octave script
    \*---------------------------------------------------------*/

    fout = fopen("tcohpsk_out.txt","wt");
    assert(fout != NULL);
    fprintf(fout, "# Created by tcohpsk.c\n");
    octave_save_int(fout, "tx_bits_log_c", tx_bits_log, 1, COHPSK_BITS_PER_FRAME*FRAMES);
    octave_save_complex(fout, "tx_symb_log_c", (COMP*)tx_symb_log, NSYMROWPILOT*FRAMES, COHPSK_NC*ND, COHPSK_NC*ND);  
    octave_save_complex(fout, "tx_fdm_frame_log_c", (COMP*)tx_fdm_frame_log, 1, COHPSK_M*NSYMROWPILOT*FRAMES, COHPSK_M*NSYMROWPILOT*FRAMES);  
    octave_save_complex(fout, "ch_fdm_frame_log_c", (COMP*)ch_fdm_frame_log_out, 1, nout-1, nout-1);  
    //octave_save_complex(fout, "rx_fdm_frame_bb_log_c", (COMP*)rx_fdm_frame_bb_log, 1, M*NSYMROWPILOT*FRAMES, M*NSYMROWPILOT*FRAMES);  
    octave_save_complex(fout, "rx_baseband_log_c", (COMP*)coh->rx_baseband_log, COHPSK_NC*ND, coh->rx_baseband_log_col_index, coh->rx_baseband_log_col_sz);  
    octave_save_complex(fout, "rx_filt_log_c", (COMP*)coh->rx_filt_log, COHPSK_NC*ND, coh->rx_filt_log_col_index, coh->rx_filt_log_col_sz);  
    octave_save_complex(fout, "ch_symb_log_c", (COMP*)coh->ch_symb_log, coh->ch_symb_log_r, COHPSK_NC*ND, COHPSK_NC*ND);  
    octave_save_float(fout, "rx_timing_log_c", (float*)coh->rx_timing_log, 1, coh->rx_timing_log_index, coh->rx_timing_log_index);  
    octave_save_complex(fout, "ct_symb_ff_log_c", (COMP*)ct_symb_ff_log, NSYMROWPILOT*FRAMES, COHPSK_NC*ND, COHPSK_NC*ND);  
    octave_save_float(fout, "rx_amp_log_c", (float*)rx_amp_log, log_data_r, COHPSK_NC*ND, COHPSK_NC*ND);  
    octave_save_float(fout, "rx_phi_log_c", (float*)rx_phi_log, log_data_r, COHPSK_NC*ND, COHPSK_NC*ND);  
    octave_save_complex(fout, "rx_symb_log_c", (COMP*)rx_symb_log, log_data_r, COHPSK_NC*ND, COHPSK_NC*ND);  
    octave_save_int(fout, "rx_bits_log_c", rx_bits_log, 1, COHPSK_BITS_PER_FRAME*log_bits);
    octave_save_float(fout, "f_est_log_c", &f_est_log[1], 1, f_est_samples-1, f_est_samples-1);  
#ifdef XX
#endif
    fclose(fout);

    cohpsk_destroy(coh);

    return 0;
}

