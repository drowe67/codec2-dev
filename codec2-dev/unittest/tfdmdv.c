/*---------------------------------------------------------------------------*\
                                                                             
  FILE........: tfdmdv.c
  AUTHOR......: David Rowe  
  DATE CREATED: April 16 2012
                                                                             
  Tests for the C version of the FDMDV modem.  This program outputs a
  file of Octave vectors that are loaded and compared to the
  Octave version of thr modem by the Octave script tfmddv.m
                                                                             
\*---------------------------------------------------------------------------*/


/*
  Copyright (C) 2009 David Rowe

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
#include "fdmdv.h"

#define FRAMES 25

void octave_save_int(FILE *f, char name[], int data[], int rows, int cols);
void octave_save_float(FILE *f, char name[], float data[], int rows, int cols);
void octave_save_complex(FILE *f, char name[], COMP data[], int rows, int cols, int col_len);

int main(int argc, char *argv[])
{
    struct FDMDV *fdmdv;
    int           tx_bits[FDMDV_BITS_PER_FRAME];
    COMP          tx_symbols[NC+1];
    COMP          tx_baseband[NC+1][M];
    COMP          tx_fdm[M];
    float         rx_fdm[M+M/P];
    float         foff;
    COMP          foff_rect;
    COMP          foff_phase_rect;
    int           nin;
    COMP          rx_fdm_fcorr[M+M/P];
    COMP          rx_baseband[NC+1][M+M/P];
    COMP          rx_filt[NC+1][P+1];
    float         rx_timing;
    float         env[NT*P];
    COMP          rx_symbols[NC+1];

    int           tx_bits_log[FDMDV_BITS_PER_FRAME*FRAMES];
    COMP          tx_symbols_log[(NC+1)*FRAMES];
    COMP          tx_baseband_log[(NC+1)][M*FRAMES];
    COMP          tx_fdm_log[M*FRAMES];
    COMP          pilot_baseband1_log[NPILOTBASEBAND*FRAMES];
    COMP          pilot_baseband2_log[NPILOTBASEBAND*FRAMES];
    COMP          pilot_lpf1_log[NPILOTLPF*FRAMES];
    COMP          pilot_lpf2_log[NPILOTLPF*FRAMES];
    COMP          S1_log[MPILOTFFT*FRAMES];
    COMP          S2_log[MPILOTFFT*FRAMES];
    float         foff_log[FRAMES];
    COMP          rx_baseband_log[(NC+1)][(M+M/P)*FRAMES];
    int           rx_baseband_log_col_index;
    COMP          rx_filt_log[NC+1][(P+1)*FRAMES];
    int           rx_filt_log_col_index;
    float         env_log[NT*P*FRAMES];
    float         rx_timing_log[FRAMES];
    COMP          rx_symbols_log[NC+1][FRAMES];
			  
    FILE         *fout;
    int           f,c,i;

    fdmdv = fdmdv_create();
    foff_phase_rect.real = 0.0; foff_phase_rect.imag = 0.0;

    rx_baseband_log_col_index = 0;
    rx_filt_log_col_index = 0;

    for(f=0; f<FRAMES; f++) {

	/* modulator -----------------------------------------*/

	fdmdv_get_test_bits(fdmdv, tx_bits);
	bits_to_dqpsk_symbols(tx_symbols, fdmdv->prev_tx_symbols, tx_bits, &fdmdv->tx_pilot_bit);
	memcpy(fdmdv->prev_tx_symbols, tx_symbols, sizeof(COMP)*(NC+1));
	tx_filter(tx_baseband, tx_symbols, fdmdv->tx_filter_memory);
	fdm_upconvert(tx_fdm, tx_baseband, fdmdv->phase_tx, fdmdv->freq);

	for(i=0; i<M; i++)
	    rx_fdm[i] = tx_fdm[i].real;
	nin = M;

	/* demodulator ----------------------------------------*/

	/* Freq offset estimation and correction */

	foff = rx_est_freq_offset(fdmdv, rx_fdm, nin);

	/* note this should be a C function with states in fdmdv */

	foff_rect.real = cos(2.0*PI*foff/FS);
	foff_rect.imag = sin(2.0*PI*foff/FS);
	for(i=0; i<nin; i++) {
	    //foff_phase_rect = cmult(foff_phase_rect, foff_rect);
	    //rx_fdm_fcorr[i] = fcmult(rx_fdm[i], foff_phase_rect);
	    rx_fdm_fcorr[i].real = rx_fdm[i];
	    rx_fdm_fcorr[i].imag = 0.0;
	}
	
	fdm_downconvert(rx_baseband, rx_fdm_fcorr, fdmdv->phase_rx, fdmdv->freq, nin);
	rx_filter(rx_filt, rx_baseband, fdmdv->rx_filter_memory, nin);
	rx_timing = rx_est_timing(rx_symbols, rx_filt, rx_baseband, fdmdv->rx_filter_mem_timing, env, fdmdv->rx_baseband_mem_timing, nin);	 

	/* save log of outputs ------------------------------------------------------*/

	memcpy(&tx_bits_log[FDMDV_BITS_PER_FRAME*f], tx_bits, sizeof(int)*FDMDV_BITS_PER_FRAME);
	memcpy(&tx_symbols_log[(NC+1)*f], tx_symbols, sizeof(COMP)*(NC+1));
	for(c=0; c<NC+1; c++)
	    for(i=0; i<M; i++)
		tx_baseband_log[c][f*M+i] = tx_baseband[c][i]; 
	memcpy(&tx_fdm_log[M*f], tx_fdm, sizeof(COMP)*M);

	/* freq offset estimation */

	memcpy(&pilot_baseband1_log[f*NPILOTBASEBAND], fdmdv->pilot_baseband1, sizeof(COMP)*NPILOTBASEBAND);
	memcpy(&pilot_baseband2_log[f*NPILOTBASEBAND], fdmdv->pilot_baseband2, sizeof(COMP)*NPILOTBASEBAND);
	memcpy(&pilot_lpf1_log[f*NPILOTLPF], fdmdv->pilot_lpf1, sizeof(COMP)*NPILOTLPF);
	memcpy(&pilot_lpf2_log[f*NPILOTLPF], fdmdv->pilot_lpf2, sizeof(COMP)*NPILOTLPF);
	memcpy(&S1_log[f*MPILOTFFT], fdmdv->S1, sizeof(COMP)*MPILOTFFT);
	memcpy(&S2_log[f*MPILOTFFT], fdmdv->S2, sizeof(COMP)*MPILOTFFT);
 	foff_log[f] = foff;

	/* rx down conversion */

	for(c=0; c<NC+1; c++) {
	    for(i=0; i<nin; i++)
		rx_baseband_log[c][rx_baseband_log_col_index + i] = rx_baseband[c][i]; 
	}
	rx_baseband_log_col_index += nin;

	/* rx filtering */

	for(c=0; c<NC+1; c++) {
	    for(i=0; i<(P*M)/nin; i++)
		rx_filt_log[c][rx_filt_log_col_index + i] = rx_filt[c][i]; 
	}
	rx_filt_log_col_index += (P*M)/nin;

	/* timing estimation */

	memcpy(&env_log[NT*P*f], env, sizeof(float)*NT*P);
	rx_timing_log[f] = rx_timing;
	for(c=0; c<NC+1; c++)
	    rx_symbols_log[c][f] = rx_symbols[c];
    }

    /* dump logs to Octave file for evaluation by tfdmdv.m Octave script */

    fout = fopen("tfdmdv_out.txt","wt");
    assert(fout != NULL);
    fprintf(fout, "# Created by tfdmdv.c\n");
    octave_save_int(fout, "tx_bits_log_c", tx_bits_log, 1, FDMDV_BITS_PER_FRAME*FRAMES);
    octave_save_complex(fout, "tx_symbols_log_c", tx_symbols_log, 1, (NC+1)*FRAMES, (NC+1)*FRAMES);  
    octave_save_complex(fout, "tx_baseband_log_c", (COMP*)tx_baseband_log, (NC+1), M*FRAMES, M*FRAMES);  
    octave_save_complex(fout, "tx_fdm_log_c", (COMP*)tx_fdm_log, 1, M*FRAMES, M*FRAMES);  
    octave_save_complex(fout, "pilot_lut_c", (COMP*)fdmdv->pilot_lut, 1, NPILOT_LUT, NPILOT_LUT);  
    octave_save_complex(fout, "pilot_baseband1_log_c", pilot_baseband1_log, 1, NPILOTBASEBAND*FRAMES, NPILOTBASEBAND*FRAMES);  
    octave_save_complex(fout, "pilot_baseband2_log_c", pilot_baseband2_log, 1, NPILOTBASEBAND*FRAMES, NPILOTBASEBAND*FRAMES);  
    octave_save_complex(fout, "pilot_lpf1_log_c", pilot_lpf1_log, 1, NPILOTLPF*FRAMES, NPILOTLPF*FRAMES);  
    octave_save_complex(fout, "pilot_lpf2_log_c", pilot_lpf2_log, 1, NPILOTLPF*FRAMES, NPILOTLPF*FRAMES);  
    octave_save_complex(fout, "S1_log_c", S1_log, 1, MPILOTFFT*FRAMES, MPILOTFFT*FRAMES);  
    octave_save_complex(fout, "S2_log_c", S2_log, 1, MPILOTFFT*FRAMES, MPILOTFFT*FRAMES);  
    octave_save_float(fout, "foff_log_c", foff_log, 1, FRAMES);  
    octave_save_complex(fout, "rx_baseband_log_c", (COMP*)rx_baseband_log, (NC+1), rx_baseband_log_col_index, (M+M/P)*FRAMES);  
    octave_save_complex(fout, "rx_filt_log_c", (COMP*)rx_filt_log, (NC+1), rx_filt_log_col_index, (P+1)*FRAMES);  
    octave_save_float(fout, "env_log_c", env_log, 1, NT*P*FRAMES);  
    octave_save_float(fout, "rx_timing_log_c", rx_timing_log, 1, FRAMES);  
    octave_save_complex(fout, "rx_symbols_log_c", (COMP*)rx_symbols_log, (NC+1), FRAMES, FRAMES);  
    fclose(fout);

    codec2_destroy(fdmdv);

    return 0;
}

void octave_save_int(FILE *f, char name[], int data[], int rows, int cols)
{
    int r,c;

    fprintf(f, "# name: %s\n", name);
    fprintf(f, "# type: matrix\n");
    fprintf(f, "# rows: %d\n", rows);
    fprintf(f, "# columns: %d\n", cols);
    
    for(r=0; r<rows; r++) {
	for(c=0; c<cols; c++)
	    fprintf(f, " %d", data[r*cols+c]);
	fprintf(f, "\n");
    }

    fprintf(f, "\n\n");
}

void octave_save_float(FILE *f, char name[], float data[], int rows, int cols)
{
    int r,c;

    fprintf(f, "# name: %s\n", name);
    fprintf(f, "# type: matrix\n");
    fprintf(f, "# rows: %d\n", rows);
    fprintf(f, "# columns: %d\n", cols);
    
    for(r=0; r<rows; r++) {
	for(c=0; c<cols; c++)
	    fprintf(f, " %f", data[r*cols+c]);
	fprintf(f, "\n");
    }

    fprintf(f, "\n\n");
}
void octave_save_complex(FILE *f, char name[], COMP data[], int rows, int cols, int col_len)
{
    int r,c;

    fprintf(f, "# name: %s\n", name);
    fprintf(f, "# type: complex matrix\n");
    fprintf(f, "# rows: %d\n", rows);
    fprintf(f, "# columns: %d\n", cols);
    printf("rows %d cols %d col_len %d\n", rows, cols, col_len);
    for(r=0; r<rows; r++) {
	for(c=0; c<cols; c++)
	    fprintf(f, " (%f,%f)", data[r*col_len+c].real, data[r*col_len+c].imag);
	fprintf(f, "\n");
    }

    fprintf(f, "\n\n");
}
