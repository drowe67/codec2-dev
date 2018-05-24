/*---------------------------------------------------------------------------*\

  FILE........: ofdm.c
  AUTHORS.....: David Rowe & Steve Sampson
  DATE CREATED: June 2017

  A Library of functions that implement a QPSK OFDM modem, C port of 
  the Octave functions in ofdm_lib.m

\*---------------------------------------------------------------------------*/

/*
  Copyright (C) 2017/2018 David Rowe

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

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <complex.h>

#include "comp.h"
#include "ofdm_internal.h"
#include "codec2_ofdm.h"
#include "ofdm_bpf_coeff.h"

/* Concrete definition of 700D parameters */
const struct OFDM_CONFIG OFDM_CONFIG_700D_C = 
{.a = 0};

/* Pointer to 700D config */ 
const struct OFDM_CONFIG  * OFDM_CONFIG_700D = &OFDM_CONFIG_700D_C;

/* Static Prototypes */

static void dft(struct OFDM *, complex float *, complex float *);
static void idft(struct OFDM *, complex float *, complex float *);
static complex float vector_sum(complex float *, int);

/* Defines */

#define max( a, b ) ( ((a) > (b)) ? (a) : (b) )
#define min( a, b ) ( ((a) < (b)) ? (a) : (b) )

/* Constants */

/*
 * QPSK Quadrant bit-pair values - Gray Coded
 */
static const complex float constellation[] = {
     1.0f + 0.0f * I,
     0.0f + 1.0f * I,
     0.0f - 1.0f * I,
    -1.0f + 0.0f * I
};

/*
 * These pilots are compatible with Octave version
 */
static const char pilotvalues[] = {
    -1, -1, 1, 1, -1, -1, -1, 1, -1,
     1, -1, 1, 1,  1,  1,  1, 1,  1,
     1
};

/*
 * A Unique Word (UW) used to verify we have modem frame sync when we
 * try a candidate coarse timing and freq offset.  The UW bits/symbols
 * are distributed through the modem frame using the index (ind)
 * tables below.  The indexes in uw_nd_sym [] and uw_ind are Octave
 * start from 1 format, subtract 1 to index C arrays.
 */

static const int tx_uw[]          = {0,0,0,0,0,0,0,0,0,0};            /* UW bits    */
static complex float tx_uw_syms[] = {1.0f,1.0f,1.0f,1.0f,1.0f};       /* UW QPSK symbols */
static const int uw_ind[]         = {19,20,37,38,55,56,73,74,91,92};  /* index in modem frame of UW bits */
static const int uw_ind_sym[]     = {10,19,28,37,46};                 /* index into modem frame of UW symbols */

/* Functions -------------------------------------------------------------------*/

/*
 * ------------
 * ofdm_create
 * ------------
 *
 * Returns OFDM data structure on success
 * Return NULL on fail
 */

struct OFDM *ofdm_create(const struct OFDM_CONFIG *config) {
    struct OFDM *ofdm;
    int i, j, n;

    if ((ofdm = (struct OFDM *) malloc(sizeof (struct OFDM))) == NULL) {
        return NULL;
    }

    /* Copy config structure */

    if (config == NULL) { /* prevent segmentation error */
        return NULL;
    }

    memcpy((void*)&ofdm->config,(void*)config,sizeof(struct OFDM_CONFIG));

    /* store complex BPSK pilot symbols */

    assert(sizeof(pilotvalues) == (OFDM_NC+2));

    for (i = 0; i < (OFDM_NC + 2); i++) {
        ofdm->pilots[i] = ((float) pilotvalues[i]) + 0.0f * I;
    }

    /* carrier tables for up and down conversion */

    float alower = OFDM_CENTRE - OFDM_RS * ((float)OFDM_NC / 2);
    int Nlower = floorf(alower / OFDM_RS);
    
    for (i = 0, n = Nlower; i < (OFDM_NC + 2); i++, n++) {
        float w = (TAU * (float) n) / (OFDM_FS / OFDM_RS);

        for (j = 0; j < OFDM_M; j++) {
            ofdm->W[i][j] = cexpf(I * w * j);
        }
    }

    for (i = 0; i < OFDM_RXBUF; i++) {
        ofdm->rxbuf[i] = 0.0f + 0.0f * I;        
    }

    for (i = 0; i < (OFDM_NS + 3); i++) {
        for (j = 0; j < (OFDM_NC + 2); j++) {
            ofdm->rx_sym[i][j] = 0.0f + 0.0f * I;
        }
    }

    for (i = 0; i < OFDM_ROWSPERFRAME*OFDM_NC; i++) {
        ofdm->rx_np[i] = 0.0f + 0.0f * I;
    }
    
    for (i = 0; i < OFDM_ROWSPERFRAME; i++) {
        for (j = 0; j < OFDM_NC; j++) {
            ofdm->aphase_est_pilot_log[OFDM_NC*i+j] = 0.0f + 0.0f * I;
            ofdm->rx_amp[OFDM_NC*i+j] = 0.0f + 0.0f * I;
        }
    }

    /* default settings of options and states */

    ofdm->verbose = 0;
    ofdm->timing_en = true;
    ofdm->foff_est_en = true;
    ofdm->phase_est_en = true;

    ofdm->foff_est_gain = 0.05f;
    ofdm->foff_est_hz = 0.0f;
    ofdm->sample_point = 0;
    ofdm->timing_est = 0;
    ofdm->timing_valid = 0;
    ofdm->timing_mx = 0.0f;
    ofdm->nin = OFDM_SAMPLESPERFRAME;
    ofdm->mean_amp = 0.0f;
    ofdm->foff_metric = 0.0f + 0.0f * I;
    
    /* sync state machine */

    for (i = 0; i < OFDM_NUWBITS; i++) {
        ofdm->tx_uw[i] = tx_uw[i];
    }

    strcpy(ofdm->sync_state,"search");
    strcpy(ofdm->last_sync_state,"search");
    ofdm->uw_errors = 0;
    ofdm->sync_counter = 0;
    ofdm->frame_count = 0;
    ofdm->sync_start = 0;
    ofdm->sync_end = 0;
    ofdm->sync_mode = OFDM_SYNC_AUTO;
    
    strcpy(ofdm->sync_state_interleaver,"search");
    strcpy(ofdm->last_sync_state_interleaver,"search");
    ofdm->frame_count_interleaver = 0;
    
    /* create the OFDM waveform */

    complex float temp[OFDM_M];

    idft(ofdm, temp, ofdm->pilots);

    /*
     * pilot_samples is 160 samples, but timing and freq offset est
     * were found by experiment to work better without a cyclic
     * prefix, so we uses zeroes instead.
     */

    /* zero out Cyclic Prefix (CP) values */

    for (i = 0, j = (OFDM_M - OFDM_NCP); i < OFDM_NCP; i++, j++) {
        ofdm->pilot_samples[i] = 0.0f + 0.0f * I;
    }
        
    // From Octave: states.timing_norm = Npsam*(rate_fs_pilot_samples*rate_fs_pilot_samples');

    /* Now copy the whole thing after the above */

    for (i = OFDM_NCP, j = 0; j < OFDM_M; i++, j++) {
        ofdm->pilot_samples[i] = temp[j];
    }

    /* calculate constant used to normalise timing correlation maximum */

    float acc = 0.0f;

    for (i = 0; i < (OFDM_M + OFDM_NCP); i++) {
        acc += crealf(ofdm->pilot_samples[i]) * crealf(ofdm->pilot_samples[i]) +
               cimagf(ofdm->pilot_samples[i]) * cimagf(ofdm->pilot_samples[i]);
    }

    ofdm->timing_norm = (OFDM_M + OFDM_NCP) * acc;

    //fprintf(stderr, "timing_norm: %f\n", ofdm->timing_norm);

    ofdm->sig_var = ofdm->noise_var = 1.0f;

    ofdm->tx_bpf_en = 0;
    ofdm->tx_bpf_buf = (complex float*)malloc(sizeof(complex float)*(OFDM_BPF_N+OFDM_SAMPLESPERFRAME));
    if (ofdm->tx_bpf_buf == NULL) {
        free(ofdm);
        return NULL;
    }
    
    for (i=0; i<OFDM_BPF_N; i++) {
        ofdm->tx_bpf_buf[i] = 0.0f + 0.0f * I;
    }
    
    return ofdm; /* Success */
}


void ofdm_destroy(struct OFDM *ofdm) {
    free(ofdm->tx_bpf_buf);
    free(ofdm);
}


/* Gray coded QPSK modulation function */

complex float qpsk_mod(int *bits) {
    return constellation[(bits[1] << 1) | bits[0]];
}

/* Gray coded QPSK demodulation function */

void qpsk_demod(complex float symbol, int *bits) {
    complex float rotate = symbol * cexpf(I * (M_PI / 4.0f));
    bits[0] = crealf(rotate) < 0.0f;
    bits[1] = cimagf(rotate) < 0.0f;
}

/* convert frequency domain into time domain */

static void idft(struct OFDM *ofdm, complex float *result, complex float *vector) {
    int row, col;

    for (row = 0; row < OFDM_M; row++) {
        result[row] = 0.0f + 0.0f * I;

        for (col = 0; col < (OFDM_NC + 2); col++) {
            result[row] = result[row] + (vector[col] * ofdm->W[col][row]);
        }

        result[row] = result[row] * OFDM_INVERSE_M;
    }
}

/* convert time domain into frequency domain */

static void dft(struct OFDM *ofdm, complex float *result, complex float *vector) {
    int row, col;

    for (col = 0; col < (OFDM_NC + 2); col++) {
        result[col] = 0.0f + 0.0f * I;

        for (row = 0; row < OFDM_M; row++) {
            result[col] = result[col] + (vector[row] * conjf(ofdm->W[col][row]));
        }
    }
}

static complex float vector_sum(complex float *a, int num_elements) {
    int i;

    complex float sum = 0.0f + 0.0f * I;

    for (i = 0; i < num_elements; i++) {
        sum = sum + a[i];
    }

    return sum;
}

/*
 * Correlates the OFDM pilot symbol samples with a window of received
 * samples to determine the most likely timing offset.  Combines two
 * frames pilots so we need at least Nsamperframe+M+Ncp samples in rx.
 *
 * Can be used for acquisition (coarse timing), and fine timing.
 *
 * Unlike Octave version use states to return a few values.
 */

static int est_timing(struct OFDM *ofdm, complex float *rx, int length) {
    complex float csam;
    int Ncorr = length - (OFDM_SAMPLESPERFRAME + (OFDM_M + OFDM_NCP));
    int SFrame = OFDM_SAMPLESPERFRAME;
    float corr[Ncorr];
    int i, j;

    float acc = 0.0f;

    for (i = 0; i < length; i++) {
        acc += crealf(rx[i]) * crealf(rx[i]) + cimagf(rx[i]) * cimagf(rx[i]);
    }
            
    float av_level = 2.0f * sqrtf(ofdm->timing_norm * acc / length) + 1E-12f;
    
    for (i = 0; i < Ncorr; i++) {
        complex float corr_st = 0.0f + 0.0f * I;
        complex float corr_en = 0.0f + 0.0f * I;

        for (j = 0; j < (OFDM_M + OFDM_NCP); j++) {
            csam = conjf(ofdm->pilot_samples[j]);
            corr_st = corr_st + (rx[i + j] * csam);
            corr_en = corr_en + (rx[i + j + SFrame] * csam);
        }

        corr[i] = (cabsf(corr_st) + cabsf(corr_en)) / av_level;
    }

    /* find the max magnitude and its index */

    float timing_mx = 0.0f;
    int timing_est = 0;

    for (i = 0; i < Ncorr; i++) {
        if (corr[i] > timing_mx) {
            timing_mx = corr[i];
            timing_est = i;
        }
    }

    ofdm->timing_mx = timing_mx;
    ofdm->timing_valid = timing_mx > OFDM_TIMING_MX_THRESH;

    if (ofdm->verbose > 1) {
        fprintf(stderr, "  av_level: %f  max: %f timing_est: %d timing_valid: %d\n", av_level, ofdm->timing_mx, timing_est, ofdm->timing_valid);
    }
    
    return timing_est;
}


/*
 * Determines frequency offset at current timing estimate, used for
 * coarse freq offset estimation during acquisition.
 *
 * Freq offset is based on an averaged statistic that was found to be
 * necessary to generate good quality estimates.
 *
 * Keep calling it when in trial or actual sync to keep statistic
 * updated, in case we lose sync.
 */

 static float est_freq_offset(struct OFDM *ofdm, complex float *rx, int length, int timing_est) {
    complex float csam1, csam2;
    int Fs = OFDM_FS;
    int SFrame = OFDM_SAMPLESPERFRAME;
    int j, k;
    float foff_est;
    
    /*
      Freq offset can be considered as change in phase over two halves
      of pilot symbols.  We average this statistic over this and next
      frames pilots.
    */

    complex float p1, p2, p3, p4;
    p1 = p2 = p3 = p4 =  0.0f + 0.0f * I;
    
    /* calculate phase of pilots at half symbol intervals */
    
    for (j = 0, k = (OFDM_M + OFDM_NCP)/2; j < (OFDM_M + OFDM_NCP)/2; j++,k++) {
        csam1 = conjf(ofdm->pilot_samples[j]);
        csam2 = conjf(ofdm->pilot_samples[k]);

        /* pilot at start of frame */
        
        p1 = p1 + (rx[timing_est + j] * csam1);
        p2 = p2 + (rx[timing_est + k] * csam2);

        /* pilot at end of frame */
        
        p3 = p3 + (rx[timing_est + j + SFrame] * csam1);
        p4 = p4 + (rx[timing_est + k + SFrame] * csam2);
    }
    
    /* Calculate sample rate of phase samples, we are sampling phase
       of pilot at half a symbol intervals */
    
    float Fs1 = Fs / ((OFDM_M + OFDM_NCP)/2);

    /* subtract phase of adjacent samples, rate of change of phase is
       frequency est.  We combine samples from either end of frame to
       improve estimate.  Small real 1E-12 term to prevent instability
       with 0 inputs. */

    ofdm->foff_metric = 0.9*ofdm->foff_metric + 0.1*(conjf(p1) * p2 + conjf(p3) * p4);
    foff_est = Fs1 * cargf( ofdm->foff_metric + 1E-12f) / TAU;

    if (ofdm->verbose > 1) {
        fprintf(stderr, "  foff_metric: %f %f foff_est: %f\n", crealf(ofdm->foff_metric), cimagf(ofdm->foff_metric), foff_est);
    }
    
    return foff_est;
}


/*
 * ----------------------------------------------
 * ofdm_txframe - modulates one frame of symbols
 * ----------------------------------------------
 */

void ofdm_txframe(struct OFDM *ofdm, complex float tx_filt[OFDM_SAMPLESPERFRAME], complex float *tx_sym_lin) {
    complex float aframe[OFDM_NS][OFDM_NC + 2];
    complex float asymbol[OFDM_M];
    complex float asymbol_cp[OFDM_M + OFDM_NCP];
    complex float tx[OFDM_SAMPLESPERFRAME];
    int i, j, k, m;

    /* initialize aframe to complex zero */

    for (i = 0; i < OFDM_NS; i++) {
        for (j = 0; j < (OFDM_NC + 2); j++) {
            aframe[i][j] = 0.0f + 0.0f * I;
        }
    }

    /* copy in a row of complex pilots to first row */

    for (i = 0; i < (OFDM_NC + 2); i++) {
        aframe[0][i] = ofdm->pilots[i];
        //printf("%2d % f % f\n", i, crealf(aframe[0][i]), cimagf(aframe[0][i]));
    }
    
    /* Place symbols in multi-carrier frame with pilots */
    /* This will place boundary values of complex zero around data */

    for (i = 1; i <= OFDM_ROWSPERFRAME; i++) {

        /* copy in the Nc complex values with [0 Nc 0] or (Nc + 2) total */

        for (j = 1; j < (OFDM_NC + 1); j++) {
            aframe[i][j] = tx_sym_lin[((i - 1) * OFDM_NC) + (j - 1)];
        }
     }

    /* OFDM up-convert symbol by symbol so we can add CP */

    for (i = 0, m = 0; i < OFDM_NS; i++, m += (OFDM_M + OFDM_NCP)) {
        idft(ofdm, asymbol, aframe[i]);

        /* Copy the last Ncp samples to the front */

        for (j = (OFDM_M - OFDM_NCP), k = 0; j < OFDM_M; j++, k++) {
            asymbol_cp[k] = asymbol[j];
        }

        /* Now copy the all samples for this row after it */

        for (j = OFDM_NCP, k = 0; k < OFDM_M; j++, k++) {
            asymbol_cp[j] = asymbol[k];
        }

        /* Now move row to the tx output */

        for (j = 0; j < (OFDM_M + OFDM_NCP); j++) {
            tx[m + j] = asymbol_cp[j];
        }
    }

    /* optional Tx Band Pass Filter */

    if (ofdm->tx_bpf_en) {
        complex float *buf = ofdm->tx_bpf_buf;
        for(i=0, j=OFDM_BPF_N; i<OFDM_SAMPLESPERFRAME; i++,j++) {
            buf[j] = tx[i];
            tx_filt[i] = 0.0;
            for(k=0; k<OFDM_BPF_N; k++) {
                tx_filt[i] += buf[j-k]*ofdm_bpf_coeff[k];
            }
        }

        assert(j <= (OFDM_BPF_N+OFDM_SAMPLESPERFRAME));
        
        /* update filter memory */

        for(i=0; i<OFDM_BPF_N; i++) {
           buf[i] = buf[i+OFDM_SAMPLESPERFRAME];
        }
    } else {
        for(i=0; i<OFDM_SAMPLESPERFRAME; i++) {
            tx_filt[i] = tx[i];
        }
    }

}

int ofdm_get_nin(struct OFDM *ofdm) {
    return ofdm->nin;
}

int ofdm_get_samples_per_frame() {
    return OFDM_SAMPLESPERFRAME;
}

int ofdm_get_max_samples_per_frame() {
    return 2*OFDM_MAX_SAMPLESPERFRAME;
}

int ofdm_get_bits_per_frame(struct OFDM *ofdm) {
    return OFDM_BITSPERFRAME;
}

void ofdm_set_verbose(struct OFDM *ofdm, int level) {
    ofdm->verbose = level;
}

void ofdm_set_timing_enable(struct OFDM *ofdm, bool val) {
    ofdm->timing_en = val;

    if (ofdm->timing_en == false) {
        /* manually set ideal timing instant */
        ofdm->sample_point = (OFDM_NCP - 1);
    }
}

void ofdm_set_foff_est_enable(struct OFDM *ofdm, bool val) {
    ofdm->foff_est_en = val;
}

void ofdm_set_phase_est_enable(struct OFDM *ofdm, bool val) {
    ofdm->phase_est_en = val;
}

void ofdm_set_off_est_hz(struct OFDM *ofdm, float val) {
    ofdm->foff_est_hz = val;
}

void ofdm_set_tx_bpf(struct OFDM *ofdm, bool val) {
    ofdm->tx_bpf_en = val;
}

/*
 * --------------------------------------
 * ofdm_mod - modulates one frame of bits
 * --------------------------------------
 */

void ofdm_mod(struct OFDM *ofdm, COMP result[OFDM_SAMPLESPERFRAME], const int *tx_bits) {
    int length = OFDM_BITSPERFRAME / OFDM_BPS;
    complex float tx[OFDM_SAMPLESPERFRAME];
    complex float tx_sym_lin[length];
    int dibit[2];
    int s, i;

    if (OFDM_BPS == 1) {
        /* Here we will have Nbitsperframe / 1 */

        for (s = 0; s < length; s++) {
            tx_sym_lin[s] = (float) (2 * tx_bits[s] - 1) + 0.0f * I;
        }
    } else if (OFDM_BPS == 2) {
        /* Here we will have Nbitsperframe / 2 */

        for (s = 0, i = 0; i < length; s += 2, i++) {
            dibit[0] = tx_bits[s + 1] & 0x1;
            dibit[1] = tx_bits[s] & 0x1;
            tx_sym_lin[i] = qpsk_mod(dibit);
        }
    }

    ofdm_txframe(ofdm, tx, tx_sym_lin);

    /* convert to comp */

    for (i = 0; i < OFDM_SAMPLESPERFRAME; i++) {
        result[i].real = crealf(tx[i]);
        result[i].imag = cimagf(tx[i]);
    }
}


/*
 * ----------------------------------------------------------------------------------
 * ofdm_sync_search - attempts to find coarse sync parameters for modem initial sync
 * ----------------------------------------------------------------------------------
 */

int ofdm_sync_search(struct OFDM *ofdm, COMP *rxbuf_in)
{
    int i,j;
    
    /* insert latest input samples into rxbuf so it is primed for when
       we have to call ofdm_demod() */

    for (i = 0, j = ofdm->nin; i < (OFDM_RXBUF - ofdm->nin); i++, j++) {
        ofdm->rxbuf[i] = ofdm->rxbuf[j];
    }

    /* insert latest input samples onto tail of rxbuf */

    for (i = (OFDM_RXBUF - ofdm->nin), j = 0; i < OFDM_RXBUF; i++, j++) {
        ofdm->rxbuf[i] = rxbuf_in[j].real + rxbuf_in[j].imag * I;
    }

    /* Attempt coarse timing estimate (i.e. detect start of frame) */

    int st = OFDM_M + OFDM_NCP + OFDM_SAMPLESPERFRAME;
    int en = st + 2*OFDM_SAMPLESPERFRAME; 
    int ct_est = est_timing(ofdm,  &ofdm->rxbuf[st], (en - st));
    ofdm->coarse_foff_est_hz = est_freq_offset(ofdm,  &ofdm->rxbuf[st], (en - st), ct_est);
   
    if (ofdm->verbose) {
        fprintf(stderr, "   ct_est: %4d foff_est: %4.1f timing_valid: %d timing_mx: %5.4f\n",
                ct_est, ofdm->coarse_foff_est_hz, ofdm->timing_valid, ofdm->timing_mx);
    }

    if (ofdm->timing_valid) {
        /* potential candidate found .... */

        /* calculate number of samples we need on next buffer to get into sync */

       ofdm->nin = OFDM_SAMPLESPERFRAME + ct_est;

       /* reset modem states */

       ofdm->sample_point = ofdm->timing_est = 0;
       ofdm->foff_est_hz = ofdm->coarse_foff_est_hz;
    } else {
        ofdm->nin = OFDM_SAMPLESPERFRAME;
    }

    return ofdm->timing_valid;
}



/*
 * ------------------------------------------
 * ofdm_demod - Demodulates one frame of bits
 * ------------------------------------------
 */

void ofdm_demod(struct OFDM *ofdm, int *rx_bits, COMP *rxbuf_in) {
    complex float aphase_est_pilot_rect;
    float aphase_est_pilot[OFDM_NC + 2];
    float aamp_est_pilot[OFDM_NC + 2];
    float freq_err_hz;
    int i, j, k, rr, st, en, ft_est;

    /* shift the buffer left based on nin */

    for (i = 0, j = ofdm->nin; i < (OFDM_RXBUF - ofdm->nin); i++, j++) {
        ofdm->rxbuf[i] = ofdm->rxbuf[j];
    }

    /* insert latest input samples onto tail of rxbuf */

    for (i = (OFDM_RXBUF - ofdm->nin), j = 0; i < OFDM_RXBUF; i++, j++) {
        ofdm->rxbuf[i] = rxbuf_in[j].real + rxbuf_in[j].imag * I;
    }

    /*
     * get user and calculated freq offset
     */

    float woff_est = TAU * ofdm->foff_est_hz / OFDM_FS;

    /* update timing estimate -------------------------------------------------- */

    if (ofdm->timing_en == true) {
        /* update timing at start of every frame */

        st = ((OFDM_M + OFDM_NCP) + OFDM_SAMPLESPERFRAME) - floorf(OFDM_FTWINDOWWIDTH / 2) + ofdm->timing_est;
        en = st + OFDM_SAMPLESPERFRAME - 1 + (OFDM_M + OFDM_NCP) + OFDM_FTWINDOWWIDTH;

        complex float work[(en - st)];

        /*
         * Adjust for the frequency error by shifting the phase
         * using a conjugate multiply
         */

        for (i = st, j = 0; i < en; i++, j++) {
            work[j] = ofdm->rxbuf[i] * cexpf(-I * woff_est * i);
        }

        ft_est = est_timing(ofdm, work, (en - st));
        ofdm->timing_est += (ft_est - ceilf(OFDM_FTWINDOWWIDTH / 2));

        /* keep the freq est statistic updated in case we lose sync,
           note we supply it with uncorrected rxbuf, note
           ofdm->coarse_fest_off_hz is unused in normal operation,
           but stored for use in tofdm.c */
    
        ofdm->coarse_foff_est_hz = est_freq_offset(ofdm, &ofdm->rxbuf[st], (en-st), ft_est);

        /* first frame in trial sync will have a better freq offset est - lets use it */

        if (ofdm->frame_count == 0) {
            ofdm->foff_est_hz = ofdm->coarse_foff_est_hz;
            woff_est = TAU * ofdm->foff_est_hz / OFDM_FS;
        }
        
        if (ofdm->verbose > 1) {
            fprintf(stderr, "  ft_est: %2d timing_est: %2d sample_point: %2d\n", ft_est, ofdm->timing_est, ofdm->sample_point);
        }

        /* Black magic to keep sample_point inside cyclic prefix.  Or something like that. */

        ofdm->sample_point = max(ofdm->timing_est + (OFDM_NCP / 4), ofdm->sample_point);
        ofdm->sample_point = min(ofdm->timing_est + OFDM_NCP, ofdm->sample_point);
    }

    /*
     * Convert the time-domain samples to the frequency-domain using the rx_sym
     * data matrix. This will be  Nc+2 carriers of 11 symbols.
     *
     * You will notice there are Nc+2 BPSK symbols for each pilot symbol, and that there
     * are Nc QPSK symbols for each data symbol.
     *
     *  XXXXXXXXXXXXXXXXX  <-- Timing Slip
     * PPPPPPPPPPPPPPPPPPP <-- Previous Frames Pilot
     *  DDDDDDDDDDDDDDDDD
     *  DDDDDDDDDDDDDDDDD
     *  DDDDDDDDDDDDDDDDD
     *  DDDDDDDDDDDDDDDDD      Ignore these past data symbols
     *  DDDDDDDDDDDDDDDDD
     *  DDDDDDDDDDDDDDDDD
     *  DDDDDDDDDDDDDDDDD
     * PPPPPPPPPPPPPPPPPPP <-- This Frames Pilot
     *  DDDDDDDDDDDDDDDDD
     *  DDDDDDDDDDDDDDDDD
     *  DDDDDDDDDDDDDDDDD
     *  DDDDDDDDDDDDDDDDD      These are the current data symbols to be decoded
     *  DDDDDDDDDDDDDDDDD
     *  DDDDDDDDDDDDDDDDD
     *  DDDDDDDDDDDDDDDDD
     * PPPPPPPPPPPPPPPPPPP <-- Next Frames Pilot
     *  DDDDDDDDDDDDDDDDD
     *  DDDDDDDDDDDDDDDDD
     *  DDDDDDDDDDDDDDDDD
     *  DDDDDDDDDDDDDDDDD      Ignore these next data symbols
     *  DDDDDDDDDDDDDDDDD
     *  DDDDDDDDDDDDDDDDD
     *  DDDDDDDDDDDDDDDDD
     * PPPPPPPPPPPPPPPPPPP <-- Future Frames Pilot
     *  XXXXXXXXXXXXXXXXX  <-- Timing Slip
     *
     * So this algorithm will have seven data symbols and four pilot symbols to process.
     * The average of the four pilot symbols is our phase estimation.
     */

    for (i = 0; i < (OFDM_NS + 3); i++) {
        for (j = 0; j < (OFDM_NC + 2); j++) {
            ofdm->rx_sym[i][j] = 0.0f + 0.0f * I;
        }
    }

    /*
     * "Previous" pilot symbol is one modem frame above.
     */

    st = (OFDM_M + OFDM_NCP) + 1 + ofdm->sample_point;
    en = st + OFDM_M;

    complex float work[OFDM_M];

    /* down-convert at current timing instant---------------------------------- */

    for (j = st, k = 0; j < en; j++, k++) {
        work[k] = ofdm->rxbuf[j] * cexpf(-I * woff_est * j);
    }

    /*
     * Each symbol is of course (OFDM_M + OFDM_NCP) samples long and
     * becomes Nc+2 carriers after DFT.
     *
     * We put this carrier pilot symbol at the top of our matrix:
     *
     * 1 .................. Nc+2
     *
     * +----------------------+
     * |    Previous Pilot    |  rx_sym[0]
     * +----------------------+
     * |                      |
     *
     */

    dft(ofdm, ofdm->rx_sym[0], work);

    /*
     * "This" pilot comes after the extra symbol alloted at the top, and after
     * the "previous" pilot and previous data symbols (let's call it, the previous
     * modem frame).
     *
     * So we will now be starting at "this" pilot symbol, and continuing to the
     * "next" pilot symbol.
     *
     * In this routine we also process the current data symbols.
     */

    for (rr = 0; rr < (OFDM_NS + 1); rr++) {
        st = (OFDM_M + OFDM_NCP) + OFDM_SAMPLESPERFRAME + (rr * (OFDM_M + OFDM_NCP)) + 1 + ofdm->sample_point;
        en = st + OFDM_M;

        /* down-convert at current timing instant---------------------------------- */

        for (j = st, k = 0; j < en; j++, k++) {
            work[k] = ofdm->rxbuf[j] * cexpf(-I * woff_est * j);
        }

        /*
         * We put these Nc+2 carrier symbols into our matrix after the previous pilot:
         *
         * 1 .................. Nc+2
         *
         * |    Previous Pilot    |  rx_sym[0]
         * +----------------------+
         * |      This Pilot      |  rx_sym[1]
         * +----------------------+
         * |         Data         |  rx_sym[2]
         * +----------------------+
         * |         Data         |  rx_sym[3]
         * +----------------------+
         * |         Data         |  rx_sym[4]
         * +----------------------+
         * |         Data         |  rx_sym[5]
         * +----------------------+
         * |         Data         |  rx_sym[6]
         * +----------------------+
         * |         Data         |  rx_sym[7]
         * +----------------------+
         * |         Data         |  rx_sym[8]
         * +----------------------+
         * |      Next Pilot      |  rx_sym[9]
         * +----------------------+
         * |                      |  rx_sym[10]
         */

        dft(ofdm, ofdm->rx_sym[rr + 1], work);
    }

    /*
     * OK, now we want to process to the "future" pilot symbol. This is after
     * the "next" modem frame.
     *
     * We are ignoring the data symbols between the "next" pilot and "future" pilot.
     * We only want the "future" pilot symbol, to perform the averaging of all pilots.
     */

    st = (OFDM_M + OFDM_NCP) + (3 * OFDM_SAMPLESPERFRAME) + 1 + ofdm->sample_point;
    en = st + OFDM_M;

    /* down-convert at current timing instant---------------------------------- */

    for (j = st, k = 0; j < en; j++, k++) {
        work[k] = ofdm->rxbuf[j] * cexpf(-I * woff_est * j);
    }

    /*
     * We put the future pilot after all the previous symbols in the matrix:
     *
     * 1 .................. Nc+2
     *
     * |                      |  rx_sym[9]
     * +----------------------+
     * |     Future Pilot     |  rx_sym[10]
     * +----------------------+
     */

    dft(ofdm, ofdm->rx_sym[OFDM_NS + 2], work);

    /*
     * We are finished now with the DFT and down conversion
     * From here on down we are in frequency domain
     */

    /* est freq err based on all carriers ------------------------------------ */

    if (ofdm->foff_est_en == true) {
        /*
         * sym[1] is (this) pilot symbol, sym[9] is (next) pilot symbol.
         *
         * By subtracting the two averages of these pilots, we find the frequency
         * by the change in phase over time.
         */

        complex float freq_err_rect = conjf(vector_sum(ofdm->rx_sym[1],
                OFDM_NC + 2)) * vector_sum(ofdm->rx_sym[OFDM_NS + 1], OFDM_NC + 2);

        /* prevent instability in atan(im/re) when real part near 0 */

        freq_err_rect = freq_err_rect + 1E-6f;

        freq_err_hz = cargf(freq_err_rect) * OFDM_RS / (TAU * OFDM_NS);
        ofdm->foff_est_hz += (ofdm->foff_est_gain * freq_err_hz);
    }

    /* OK - now estimate and correct pilot phase  ---------------------------------- */

    for (i = 0; i < (OFDM_NC + 2); i++) {
        aphase_est_pilot[i] = 10.0f;
        aamp_est_pilot[i] = 0.0f;
    }

    /*
     * Basically we divide the Nc+2 pilots into groups of 3
     *
     * Then average the phase surrounding each of the data symbols.
     */

    for (i = 1; i < (OFDM_NC + 1); i++) {
        complex float symbol[3];

        for (j = (i - 1), k = 0; j < (i + 2); j++, k++) {
            symbol[k] = ofdm->rx_sym[1][j] * conjf(ofdm->pilots[j]); /* this pilot conjugate */
        }

        aphase_est_pilot_rect = vector_sum(symbol, 3);

        for (j = (i - 1), k = 0; j < (i + 2); j++, k++) {
            symbol[k] = ofdm->rx_sym[OFDM_NS + 1][j] * conjf(ofdm->pilots[j]); /* next pilot conjugate */
        }

        aphase_est_pilot_rect = aphase_est_pilot_rect + vector_sum(symbol, 3);

        /* use next step of pilots in past and future */

        for (j = (i - 1), k = 0; j < (i + 2); j++, k++) {
            symbol[k] = ofdm->rx_sym[0][j] * ofdm->pilots[j]; /* previous pilot */
        }

        aphase_est_pilot_rect = aphase_est_pilot_rect + vector_sum(symbol, 3);

        for (j = (i - 1), k = 0; j < (i + 2); j++, k++) {
            symbol[k] = ofdm->rx_sym[OFDM_NS + 2][j] * ofdm->pilots[j]; /* last pilot */
        }

        aphase_est_pilot_rect = aphase_est_pilot_rect + vector_sum(symbol, 3);
        aphase_est_pilot[i] = cargf(aphase_est_pilot_rect);

        /* amplitude is estimated over 12 pilots */

        aamp_est_pilot[i] = cabsf(aphase_est_pilot_rect / 12.0f);
    }

    /*
     * correct phase offset using phase estimate, and demodulate
     * bits, separate loop as it runs across cols (carriers) to get
     * frame bit ordering correct
     */

    complex float rx_corr;
    int abit[2];
    int bit_index = 0;
    float sum_amp = 0.0f;
    
    for (rr = 0; rr < OFDM_ROWSPERFRAME; rr++) {
        /*
         * Note the i starts with the second carrier, ends with Nc+1.
         * so we ignore the first and last carriers.
         *
         * Also note we are using sym[2..8] or the seven data symbols.
         */

        for (i = 1; i < (OFDM_NC + 1); i++) {
            if (ofdm->phase_est_en == true) {
                rx_corr = ofdm->rx_sym[rr + 2][i] * cexpf(-I * aphase_est_pilot[i]);
            } else {
                rx_corr = ofdm->rx_sym[rr + 2][i];
            }

            /*
             * Output complex data symbols after phase correction;
             * rx_np means the pilot symbols have been removed
             */

            ofdm->rx_np[(rr * OFDM_NC) + (i - 1)] = rx_corr;

            /*
             * Note even though amp ests are the same for each col,
             * the FEC decoder likes to have one amplitude per symbol
             * so convenient to log them all
             */

            ofdm->rx_amp[(rr * OFDM_NC) + (i - 1)] = aamp_est_pilot[i];
            sum_amp += aamp_est_pilot[i];
            
            /*
             * Note like amps in this implementation phase ests are the
             * same for each col, but we log them for each symbol anyway
             */

            ofdm->aphase_est_pilot_log[(rr * OFDM_NC) + (i - 1)] = aphase_est_pilot[i];

            if (OFDM_BPS == 1) {
                rx_bits[bit_index++] = crealf(rx_corr) > 0.0f;
            } else if (OFDM_BPS == 2) {
                /*
                 * Only one final task, decode what quadrant the phase
                 * is in, and return the dibits
                 */
                qpsk_demod(rx_corr, abit);
                rx_bits[bit_index++] = abit[1];
                rx_bits[bit_index++] = abit[0];
            }
        }
    }

    /* update mean amplitude estimate for LDPC decoder scaling */
    
    ofdm->mean_amp = 0.9f * ofdm->mean_amp + 0.1f * sum_amp/(OFDM_ROWSPERFRAME * OFDM_NC);
    
    /* Adjust nin to take care of sample clock offset */

    ofdm->nin = OFDM_SAMPLESPERFRAME;

    if (ofdm->timing_en == true) {
        int thresh = (OFDM_M + OFDM_NCP) / 8;
        int tshift = (OFDM_M + OFDM_NCP) / 4;

        if (ofdm->timing_est > thresh) {
            ofdm->nin = OFDM_SAMPLESPERFRAME + tshift;
            ofdm->timing_est -= tshift;
            ofdm->sample_point -= tshift;
        } else if (ofdm->timing_est < -thresh) {
            ofdm->nin = OFDM_SAMPLESPERFRAME - tshift;
            ofdm->timing_est += tshift;
            ofdm->sample_point += tshift;
        }
    }

    /* estimate signal and noise power, see ofdm_lib.m, cohpsk.m for
       more info */

    float sig_var = 0.0f;
    complex float *rx_np = ofdm->rx_np;

    for (i = 0; i < (OFDM_ROWSPERFRAME * OFDM_NC); i++) {
        sig_var += crealf(rx_np[i]) * crealf(rx_np[i]) + cimagf(rx_np[i]) * cimagf(rx_np[i]);
    }

    sig_var /= (OFDM_ROWSPERFRAME * OFDM_NC);
    float sig_rms = sqrtf(sig_var);

    complex float s;
    float sum_x = 0.0f;
    float sum_xx = 0.0f;
    int n = 0;

    for (i=0; i<OFDM_ROWSPERFRAME * OFDM_NC; i++) {
      s = rx_np[i];

      if (fabsf(crealf(s)) > sig_rms) {
          sum_x  += cimagf(s);
          sum_xx += cimagf(s) * cimagf(s);
          n++;
      }
    }

    /* with large interfering carriers this alg can break down - in
       that case set a benign value for noise_var that will produce a
       sensible (probably low) SNR est */
    
    float noise_var = 1.0f;

    if (n > 1) {
        noise_var = (n * sum_xx - sum_x * sum_x) / (n * (n - 1));
    }

    ofdm->noise_var = 2.0f * noise_var;
    ofdm->sig_var = sig_var;
}


/* iterate state machine ------------------------------------*/

void ofdm_sync_state_machine(struct OFDM *ofdm, int *rx_uw) {
    char next_state[OFDM_STATE_STR];
    int  i;
    
    strcpy(next_state, ofdm->sync_state);    
    ofdm->sync_start = ofdm->sync_end = 0;
  
    if (strcmp(ofdm->sync_state,"search") == 0) { 
        if (ofdm->timing_valid) {
            ofdm->frame_count = 0;
            ofdm->sync_counter = 0;
            ofdm->sync_start = 1;
            strcpy(next_state, "trial");
        }
    }

    if (!strcmp(ofdm->sync_state,"synced") || !strcmp(ofdm->sync_state, "trial")) {
        ofdm->frame_count++;
        ofdm->frame_count_interleaver++;
       
        /* freq offset est may be too far out, and has aliases every 1/Ts, so
           we use a Unique Word to get a really solid indication of sync. */

        ofdm->uw_errors = 0;

        for (i=0; i<OFDM_NUWBITS; i++) {
            ofdm->uw_errors += ofdm->tx_uw[i] ^ rx_uw[i]; 
        }

        /* during trial sync we don't tolerate errors so much, we look
           for 3 consecutive frames with low error rate to confirm
           sync */
      
        if (!strcmp(ofdm->sync_state, "trial")) {
            if (ofdm->uw_errors > 2) {
                /* if we exceed thresh stay in trial sync */
                ofdm->sync_counter++;
                ofdm->frame_count = 0; 
            }

            if (ofdm->sync_counter == 2) {
                /* if we get two bad frames drop sync and start again */
                strcpy(next_state, "search");
                strcpy(ofdm->sync_state_interleaver, "search");                
            }
           
            if (ofdm->frame_count == 4) {
                /* three good frames, sync is OK! */
                strcpy(next_state, "synced");
            }
        }

        /* once we have synced up we tolerate a higher error rate to wait out fades */

        if (!strcmp(ofdm->sync_state, "synced")) {
            if (ofdm->uw_errors > 2) {
                ofdm->sync_counter++;
            } else {
                ofdm->sync_counter = 0;
            }
                
            if ((ofdm->sync_mode == OFDM_SYNC_AUTO) && (ofdm->sync_counter == 12)) {
                /* run of consective bad frames ... drop sync */
                strcpy(next_state, "search");
                strcpy(ofdm->sync_state_interleaver, "search");
            }           
        }
    }
    
    strcpy(ofdm->last_sync_state, ofdm->sync_state);
    strcpy(ofdm->last_sync_state_interleaver, ofdm->sync_state_interleaver);
    strcpy(ofdm->sync_state, next_state);
}


/*---------------------------------------------------------------------------* \

  FUNCTIONS...: ofdm_set_sync
  AUTHOR......: David Rowe
  DATE CREATED: May 2018

  Operator control of sync state machine.  This mode is required to
  acquire sync up at very low SNRS.  This is difficult to implement,
  for example we may get a false sync, or the state machine may fall
  out of sync by mistake during a long fade.

  So with this API call we allow some operator assistance.

  Ensure this is called in the same thread as ofdm_sync_state_machine().

\*---------------------------------------------------------------------------*/

void ofdm_set_sync(struct OFDM *ofdm, int sync_cmd) {
    assert (ofdm != NULL);

    switch(sync_cmd) {
    case OFDM_SYNC_UNSYNC:
        /* force manual unsync, in case operator detects false sync,
           which will cuase sync state machine to have another go at
           sync */
        strcpy(ofdm->sync_state, "search");
        strcpy(ofdm->sync_state_interleaver, "search");                
        break;
    case OFDM_SYNC_AUTO:
        /* normal operating mode - sync state machine decides when to unsync */
        ofdm->sync_mode = OFDM_SYNC_AUTO;
        break;
    case OFDM_SYNC_MANUAL:
        /* allow sync state machine to sync, but not to unsync, the
           operator will decide that manually */
        ofdm->sync_mode = OFDM_SYNC_MANUAL;
        break;
    default:
        assert(0);
    }            
}

/*---------------------------------------------------------------------------*\

  FUNCTION....: ofdm_get_demod_stats()
  AUTHOR......: David Rowe
  DATE CREATED: May 2018

  Fills stats structure with a bunch of demod information.

\*---------------------------------------------------------------------------*/

void ofdm_get_demod_stats(struct OFDM *ofdm, struct MODEM_STATS *stats)
{
    int   c, r;

    stats->Nc = OFDM_NC;
    assert(stats->Nc <= MODEM_STATS_NC_MAX);

    float snr_est = 10.0f * log10f((0.1+ (ofdm->sig_var/ofdm->noise_var)) * OFDM_NC*OFDM_RS / 3000.0f);
    //fprintf(stderr, "sig: %f var: %f snr: %f\n", ofdm->sig_var, ofdm->noise_var, snr_est);
    stats->snr_est = 0.9f * stats->snr_est + 0.1f * snr_est;
    stats->sync = !strcmp(ofdm->sync_state, "synced") || !strcmp(ofdm->sync_state, "trial");
    //fprintf(stderr, "sync: %d %s\n", stats->sync, ofdm->sync_state);
    stats->foff = ofdm->foff_est_hz;
    stats->rx_timing = ofdm->timing_est;
    stats->clock_offset = 0.0f;               /* TODO: work out sample clock offset */

    assert(OFDM_ROWSPERFRAME < MODEM_STATS_NR_MAX);
    stats->nr = OFDM_ROWSPERFRAME;

    for (c = 0; c < OFDM_NC; c++) {
        for (r = 0; r < OFDM_ROWSPERFRAME; r++) {
            complex float rot = ofdm->rx_np[r * c] * cexpf(I * (M_PI / 4.0f));
            stats->rx_symbols[r][c].real = crealf(rot);
            stats->rx_symbols[r][c].imag = cimagf(rot);
            //fprintf(stderr, "%f %f\n", stats->rx_symbols[r][c].real, stats->rx_symbols[r][c].imag);
        }
    }
}


/* Assemble modem frame from UW, payload symbols, and txt bits */

void ofdm_assemble_modem_frame(complex float modem_frame[],
                               COMP    payload_syms[],
                               uint8_t txt_bits[])
{
  int Nsymsperframe = OFDM_BITSPERFRAME/OFDM_BPS;
  int Nuwsyms = OFDM_NUWBITS/OFDM_BPS;
  int Ntxtsyms = OFDM_NTXTBITS/OFDM_BPS;

  int s,p=0,u=0;

  for (s=0; s<Nsymsperframe-Ntxtsyms; s++) {
      if ((u < Nuwsyms) && (s == (uw_ind_sym[u]-1))) {
          modem_frame[s] = tx_uw_syms[u++];
      } else {
          modem_frame[s] = payload_syms[p].real + payload_syms[p].imag * I;
          p++;
      }
  }
  assert(u == Nuwsyms);
  assert(p == (Nsymsperframe-Nuwsyms-Ntxtsyms));

  int t; int dibit[2];

  for (t=0; s<Nsymsperframe; s++,t+=OFDM_BPS) {
      dibit[0] = txt_bits[t+1] & 0x1;
      dibit[1] = txt_bits[t] & 0x1;
      modem_frame[s] = qpsk_mod(dibit);
  }
  assert(t == OFDM_NTXTBITS);
}


void ofdm_disassemble_modem_frame(struct OFDM   *ofdm,
                                  int            rx_uw[],
                                  COMP           codeword_syms[],
                                  float          codeword_amps[],
                                  short          txt_bits[])
{
  int Nsymsperframe = OFDM_BITSPERFRAME/OFDM_BPS;
  int Nuwsyms = OFDM_NUWBITS/OFDM_BPS;
  int Ntxtsyms = OFDM_NTXTBITS/OFDM_BPS;
  int dibit[2];
  
  int s,p=0,u=0;

  for (s=0; s<Nsymsperframe-Ntxtsyms; s++) {
      if ((u < Nuwsyms) && (s == (uw_ind_sym[u]-1))) {
          qpsk_demod(ofdm->rx_np[s], dibit);
          rx_uw[OFDM_BPS*u]   = dibit[1];
          rx_uw[OFDM_BPS*u+1] = dibit[0];
          u++;
      } else {
          codeword_syms[p].real = crealf(ofdm->rx_np[s]);
          codeword_syms[p].imag = cimagf(ofdm->rx_np[s]);
          codeword_amps[p] = ofdm->rx_amp[s];
          p++;
      }
  }
  assert(u == Nuwsyms);
  assert(p == (Nsymsperframe-Nuwsyms-Ntxtsyms));

  int t;

  for (t=0; s<Nsymsperframe; s++,t+=OFDM_BPS) {
      qpsk_demod(ofdm->rx_np[s], dibit);
      txt_bits[t]   = dibit[1];
      txt_bits[t+1] = dibit[0];
  }
  assert(t == OFDM_NTXTBITS);
}

