/*---------------------------------------------------------------------------*\

  FILE........: ofdm.c
  AUTHORS.....: David Rowe & Steve Sampson
  DATE CREATED: June 2017

  A Library of functions that implement a BPSK/QPSK OFDM modem

\*---------------------------------------------------------------------------*/

/*
  Copyright (C) 2017 David Rowe

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
#include "kiss_fft.h"

/* Concrete definition of 700D parameters */
const struct OFDM_CONFIG OFDM_CONFIG_700D_C = 
{.a = 0};

/* Pointer to 700D config */ 
const struct OFDM_CONFIG  * OFDM_CONFIG_700D = &OFDM_CONFIG_700D_C;

/* Static Prototypes */

static void dft(struct OFDM *, complex float *, complex float *);
static void idft(struct OFDM *, complex float *, complex float *);
static complex float vector_sum(complex float *, int);
static complex float qpsk_mod(int *);
static void qpsk_demod(complex float, int *);
static void ofdm_txframe(struct OFDM *, complex float [OFDM_SAMPLESPERFRAME], complex float *);
static int coarse_sync(struct OFDM *, complex float *, int, float *foff_est);

/* Defines */

#define max( a, b ) ( ((a) > (b)) ? (a) : (b) )
#define min( a, b ) ( ((a) < (b)) ? (a) : (b) )

/* Constants */

/*
 * QPSK Quadrant bit-pair values - Gray Coded
 *
 *   0.0 -  89.9 = 00
 *  90.0 - 179.9 = 01
 * 180.0 - 269.9 = 11
 * 270.0 - 359.9 = 10
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
     1, -1, 1, 1,  1,  1,  1, 1,  1
};

/* Functions */

/* Gray coded QPSK modulation function */

static complex float qpsk_mod(int *bits) {
    return constellation[(bits[1] << 1) | bits[0]];
}

/* Gray coded QPSK demodulation function */

static void qpsk_demod(complex float symbol, int *bits) {
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
            result[row] = result[row] + (vector[col] * (ofdm->W[col][row] / (float) OFDM_M)); /* complex result */
        }
    }
}

/* convert time domain into frequency domain */

static void dft(struct OFDM *ofdm, complex float *result, complex float *vector) {
    int row, col;

    for (col = 0; col < (OFDM_NC + 2); col++) {
        result[col] = 0.0f + 0.0f * I;

        for (row = 0; row < OFDM_M; row++) {
            result[col] = result[col] + (vector[row] * conjf(ofdm->W[col][row])); /* complex result */
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
 * Also estimates coarse frequency offset based on sampling the pilots
 * phase at half symbol intervales.
 * 
 * Unlike Octave version use states to return a few values.
 */

static int coarse_sync(struct OFDM *ofdm, complex float *rx, int length, float *foff_est) {
    complex float csam, csam1, csam2;
    int Ncorr = length - (OFDM_SAMPLESPERFRAME + (OFDM_M + OFDM_NCP));
    int Fs = OFDM_FS;
    int SFrame = OFDM_SAMPLESPERFRAME;
    float corr[Ncorr], av_level;
    int i, j, k;

    complex float acc =  0.0f + 0.0f * I;
    for (i = 0; i < length; i++) {
        acc += crealf(rx[i]) * crealf(rx[i]) + cimagf(rx[i]) * cimagf(rx[i]);
    }
            
    av_level = 2.0*sqrt(ofdm->timing_norm*crealf(acc)/length) + 1E-12;
    
    for (i = 0; i < Ncorr; i++) {
        complex float temp = 0.0f + 0.0f * I;

        for (j = 0; j < (OFDM_M + OFDM_NCP); j++) {
            csam = conjf(ofdm->pilot_samples[j]);
            temp = temp + (rx[i + j] * csam);
            temp = temp + (rx[i + j + SFrame] * csam);
        }

        corr[i] = cabsf(temp)/av_level;
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
        fprintf(stderr, "   max: %f timing_est: %d timing_valid: %d\n", ofdm->timing_mx, timing_est, ofdm->timing_valid);
    }
    
    /* Coarse frequency estimation - note this isn't always used, some CPU could be saved
       buy putting this in another function */

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
    
    float Fs1 = Fs/((OFDM_M + OFDM_NCP)/2);

    /* subtract phase of adjacent samples, rate of change of phase is
       frequency est.  We combine samples from either end of frame to
       improve estimate.  Small real 1E-12 term to prevent instability
       with 0 inputs. */
    
    *foff_est = Fs1 * cargf(conjf(p1)*p2 + conjf(p3)*p4 + 1E-12)/(2.0*M_PI);

    return timing_est;
}

/*
 * ----------------------------------------------
 * ofdm_txframe - modulates one frame of symbols
 * ----------------------------------------------
 */

static void ofdm_txframe(struct OFDM *ofdm, complex float tx[OFDM_SAMPLESPERFRAME],
    complex float *tx_sym_lin) {
    complex float aframe[OFDM_NS][OFDM_NC + 2];
    complex float asymbol[OFDM_M];
    complex float asymbol_cp[OFDM_M + OFDM_NCP];
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

        /* Copy the last Ncp columns to the front */

        for (j = (OFDM_M - OFDM_NCP), k = 0; j < OFDM_M; j++, k++) {
            asymbol_cp[k] = asymbol[j];
        }

        /* Now copy the whole row after it */

        for (j = OFDM_NCP, k = 0; k < OFDM_M; j++, k++) {
            asymbol_cp[j] = asymbol[k];
        }

        /* Now move row to the tx reference */

        for (j = 0; j < (OFDM_M + OFDM_NCP); j++) {
            tx[m + j] = asymbol_cp[j];
        }
    }
}

/*
 * ------------
 * ofdm_create
 * ------------
 *
 * Returns OFDM data structure on success
 * Return NULL on fail
 */

struct OFDM *ofdm_create(const struct OFDM_CONFIG * config) {
    struct OFDM *ofdm;
    int i, j;

    if ((ofdm = (struct OFDM *) malloc(sizeof (struct OFDM))) == NULL) {
        return NULL;
    }

    /* Copy config structure */
    memcpy((void*)&ofdm->config,(void*)config,sizeof(struct OFDM_CONFIG));

    /* store complex BPSK pilot symbols */

    for (i = 0; i < (OFDM_NC + 2); i++) {
        ofdm->pilots[i] = ((float) pilotvalues[i]) + 0.0f * I;
    }

    /* carrier tables for up and down conversion */

    int Nlower = floorf((OFDM_CENTRE - OFDM_RS * (OFDM_NC / 2)) / OFDM_RS);

    for (i = 0, j = Nlower; i < (OFDM_NC + 2); i++, j++) {
        /*
         * 2 * pi * j/144  j=19..36
         * j = 1 kHz to 2 kHz (1.5 kHz center)
         */

        ofdm->w[i] = TAU * (float) j / (OFDM_FS / OFDM_RS);
    }

    for (i = 0; i < (OFDM_NC + 2); i++) {
        for (j = 0; j < OFDM_M; j++) {
            ofdm->W[i][j] = cexpf(I * ofdm->w[i] * j);
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

    ofdm->foff_est_gain = 0.01f;
    ofdm->foff_est_hz = 0.0f;
    ofdm->sample_point = 0;
    ofdm->timing_est = 0;
    ofdm->timing_valid = 0;
    ofdm->timing_mx = 0.0f;
    ofdm->nin = OFDM_SAMPLESPERFRAME;

    /* create the OFDM waveform */

    complex float temp[OFDM_M];

    idft(ofdm, temp, ofdm->pilots);

    /*
     * pilot_samples is 160 samples, as we copy the last 16 to the front
     */

    /* first copy the last Cyclic Prefix (CP) values */

    for (i = 0, j = (OFDM_M - OFDM_NCP); i < OFDM_NCP; i++, j++) {
        ofdm->pilot_samples[i] = temp[j];
    }
        
    // states.timing_norm = Npsam*(rate_fs_pilot_samples*rate_fs_pilot_samples');

    /* Now copy the whole thing after the above */

    for (i = OFDM_NCP, j = 0; j < OFDM_M; i++, j++) {
        ofdm->pilot_samples[i] = temp[j];
    }

    /* calculate constant used to normalise timing correlation maximum */

    complex float acc =  0.0f + 0.0f * I;
    for (i = 0; i < OFDM_M+OFDM_NCP; i++) {
        acc += ofdm->pilot_samples[i] * conjf(ofdm->pilot_samples[i]);
    }
    ofdm->timing_norm = (OFDM_M + OFDM_NCP) * crealf(acc);
    //fprintf(stderr, "timing_norm: %f\n", ofdm->timing_norm);
    
    return ofdm; /* Success */
}

void ofdm_destroy(struct OFDM *ofdm) {
    free(ofdm);
}

int ofdm_get_nin(struct OFDM *ofdm) {
    return ofdm->nin;
}

int ofdm_get_samples_per_frame() {
    return OFDM_SAMPLESPERFRAME;
}

int ofdm_get_max_samples_per_frame() {
    return OFDM_MAX_SAMPLESPERFRAME;
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
    int ct_est = coarse_sync(ofdm,  &ofdm->rxbuf[st], (en - st), &ofdm->coarse_foff_est_hz);
    if (ofdm->verbose) {
        fprintf(stderr, "   ct_est: %4d foff_est: %3.1f timing_valid: %d timing_mx: %f\n",
                ct_est, ofdm->coarse_foff_est_hz, ofdm->timing_valid, ofdm->timing_mx);
    }

    if (ofdm->timing_valid) {
        /* potential candidate found .... */

        /* calculate number of samples we need on next buffer to get into sync */

       ofdm->nin = OFDM_SAMPLESPERFRAME + ct_est;

       /* reset modem states */

       ofdm->sample_point = ofdm->timing_est = 0;
       ofdm->foff_est_hz = ofdm->coarse_foff_est_hz;
    }
    else {
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
    float freq_err_hz, coarse_foff_est;
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

        /* note coarse sync just used for timing est, we dont use coarse_foff_est in this call */
        
        ft_est = coarse_sync(ofdm, work, (en - st), &ofdm->coarse_foff_est_hz);
        ofdm->timing_est += (ft_est - ceilf(OFDM_FTWINDOWWIDTH / 2));

        if (ofdm->verbose > 1) {
            fprintf(stdout, "  ft_est: %2d timing_est: %2d sample_point: %2d\n", ft_est, ofdm->timing_est, ofdm->sample_point);
        }

        /* Black magic to keep sample_point inside cyclic prefix.  Or something like that. */

        ofdm->sample_point = max(ofdm->timing_est + (OFDM_NCP / 4), ofdm->sample_point);
        ofdm->sample_point = min(ofdm->timing_est + OFDM_NCP, ofdm->sample_point);
    }

    /*
     * Convert the time-domain samples to the frequency-domain using the rx_sym
     * data matrix. This will be 18 carriers of 11 symbols.
     *
     * You will notice there are 18 BPSK symbols for each pilot symbol, and that there
     * are 16 QPSK symbols for each data symbol.
     *
     *  XXXXXXXXXXXXXXXX  <-- Timing Slip
     * PPPPPPPPPPPPPPPPPP <-- Previous Frames Pilot
     *  DDDDDDDDDDDDDDDD
     *  DDDDDDDDDDDDDDDD
     *  DDDDDDDDDDDDDDDD
     *  DDDDDDDDDDDDDDDD      Ignore these past data symbols
     *  DDDDDDDDDDDDDDDD
     *  DDDDDDDDDDDDDDDD
     *  DDDDDDDDDDDDDDDD
     * PPPPPPPPPPPPPPPPPP <-- This Frames Pilot
     *  DDDDDDDDDDDDDDDD
     *  DDDDDDDDDDDDDDDD
     *  DDDDDDDDDDDDDDDD
     *  DDDDDDDDDDDDDDDD      These are the current data symbols to be decoded
     *  DDDDDDDDDDDDDDDD
     *  DDDDDDDDDDDDDDDD
     *  DDDDDDDDDDDDDDDD
     * PPPPPPPPPPPPPPPPPP <-- Next Frames Pilot
     *  DDDDDDDDDDDDDDDD
     *  DDDDDDDDDDDDDDDD
     *  DDDDDDDDDDDDDDDD
     *  DDDDDDDDDDDDDDDD      Ignore these next data symbols
     *  DDDDDDDDDDDDDDDD
     *  DDDDDDDDDDDDDDDD
     *  DDDDDDDDDDDDDDDD
     * PPPPPPPPPPPPPPPPPP <-- Future Frames Pilot
     *  XXXXXXXXXXXXXXXX  <-- Timing Slip
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
     * becomes 18 carriers after DFT.
     *
     * We put this 18 carrier pilot symbol at the top of our matrix:
     *
     * 0 ................... 17
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
         * We put these 18 carrier symbols into our matrix after the previous pilot:
         *
         * 0 ................... 17
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
     * 0 ................... 17
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
     * Basically we divide the 18 pilots into groups of 3
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

        /* TODO David: WTF 12.0 constant?  Something to do with LDPC input scaling? */

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

    for (rr = 0; rr < OFDM_ROWSPERFRAME; rr++) {
        /*
         * Note the i has an index of 1 to 16, so we ignore carriers 0 and 17.
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
             * For testing, we want to save these complex data values
             * without the pilots. Thus, the name rx (no pilot) np
             */

            ofdm->rx_np[(rr * OFDM_NC) + (i - 1)] = rx_corr;

            /*
             * Note even though amp ests are the same for each col,
             * the FEC decoder likes to have one amplitude per symbol
             * so convenient to log them all
             */

            ofdm->rx_amp[(rr * OFDM_NC) + (i - 1)] = aamp_est_pilot[i];

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
}

