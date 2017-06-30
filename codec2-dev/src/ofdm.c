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

/* Static Prototypes */

static void matrix_vector_multiply(struct OFDM *, complex float *, complex float *);
static void matrix_vector_conjugate_multiply(struct OFDM *, complex float *, complex float *);
static complex float vector_sum(complex float *, int);
static complex float qpsk_mod(int *);
static void qpsk_demod(complex float, int *);
static void ofdm_txframe(struct OFDM *, complex float [OFDM_SAMPLESPERFRAME], complex float *);
static int coarse_sync(struct OFDM *, complex float *, int);

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
    1, -1, 1, 1, 1, 1, 1, 1, 1
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

static void matrix_vector_multiply(struct OFDM *ofdm, complex float *result, complex float *vector) {
    int row, col;

    for (row = 0; row < OFDM_M; row++) {
        result[row] = 0.0f + 0.0f * I;

        for (col = 0; col < (OFDM_NC + 2); col++) {
            result[row] = result[row] + (vector[col] * (ofdm->W[col][row] / (float) OFDM_M)); /* complex result */
        }
    }
}

/* convert time domain into frequency domain */

static void matrix_vector_conjugate_multiply(struct OFDM *ofdm, complex float *result, complex float *vector) {
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
 */

static int coarse_sync(struct OFDM *ofdm, complex float *rx, int length) {
    complex float csam;
    int Ncorr = length - (OFDM_SAMPLESPERFRAME + (OFDM_M + OFDM_NCP));
    float corr[Ncorr];
    int i, j;

    for (i = 0; i < Ncorr; i++) {
        complex float temp = 0.0f + 0.0f * I;

        for (j = 0; j < (OFDM_M + OFDM_NCP); j++) {
            csam = conjf(ofdm->rate_fs_pilot_samples[j]);
            temp = temp + (rx[i + j] * csam);
            temp = temp + (rx[i + j + OFDM_SAMPLESPERFRAME] * csam);
        }

        corr[i] = cabsf(temp);
    }

    /* find the max magnitude and its index */

    float mag = 0.0f;
    int t_est = 0;

    for (i = 0; i < Ncorr; i++) {
        if (corr[i] > mag) {
            mag = corr[i];
            t_est = i;
        }
    }

    return t_est;
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
        matrix_vector_multiply(ofdm, asymbol, aframe[i]);

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

struct OFDM *ofdm_create() {
    struct OFDM *ofdm;
    int i, j;

    if ((ofdm = (struct OFDM *) malloc(sizeof (struct OFDM))) == NULL) {
        return NULL;
    }

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

    for (i = 0; i < (OFDM_NS + 3); i++) {
        for (j = 0; j < (OFDM_NC + 2); j++) {
            ofdm->rx_sym[i][j] = 0.0f + 0.0f * I;
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
    ofdm->nin = OFDM_SAMPLESPERFRAME;

    /* create the OFDM waveform */

    complex float temp[OFDM_M];

    matrix_vector_multiply(ofdm, temp, ofdm->pilots);

    /*
     * rate_fs_pilot_samples is 160 samples, as we copy the last 16 to the front
     */

    /* first copy the last Cyclic Prefix (CP) values */

    for (i = 0, j = (OFDM_M - OFDM_NCP); i < OFDM_NCP; i++, j++) {
        ofdm->rate_fs_pilot_samples[i] = temp[j];
    }

    /* Now copy the whole thing after the above */

    for (i = OFDM_NCP, j = 0; j < OFDM_M; i++, j++) {
        ofdm->rate_fs_pilot_samples[i] = temp[j];
    }

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
 * ------------------------------------------
 * ofdm_demod - Demodulates one frame of bits
 * ------------------------------------------
 *
 * For phase estimation we need to maintain buffer of 3 modem frames
 * plus one pilot, so we have 4 pilots total. '^' is the start of
 * current frame that we are demodulating.
 *
 * P DDDDDDD P DDDDDDD P DDDDDDD P
 *           ^
 *
 * Then add one data symbol (M+Ncp) either side to account for
 * movement in sampling instant due to sample clock differences:
 *
 * D P DDDDDDD P DDDDDDD P DDDDDDD P D
 *             ^
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
     * get latest freq offset estimate
     *
     * foff_est_hz will be 0.0 unless
     * foff_est_en is enabled
     */

    float woff_est = TAU * ofdm->foff_est_hz / OFDM_FS;

    /* update timing estimate -------------------------------------------------- */

    if (ofdm->timing_en == true) {
        /* update timing at start of every frame */

        st = (OFDM_M + OFDM_NCP + OFDM_SAMPLESPERFRAME) - floorf(OFDM_FTWINDOWWIDTH / 2) + ofdm->timing_est;
        en = st + OFDM_SAMPLESPERFRAME - 1 + OFDM_M + OFDM_NCP + OFDM_FTWINDOWWIDTH;

        complex float work[(en - st)];

        for (i = st, j = 0; i < en; i++, j++) {
            work[j] = ofdm->rxbuf[i] * cexpf(-I * woff_est * i);
        }

        ft_est = coarse_sync(ofdm, work, (en - st));
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
     * data matrix. This will be 18 carriers of 11 symbols (P P DDDDDDD P P)
     *
     * So we will have one modem data frame and three pilots to do magic.
     */

    for (i = 0; i < (OFDM_NS + 3); i++) {
        for (j = 0; j < (OFDM_NC + 2); j++) {
            ofdm->rx_sym[i][j] = 0.0f + 0.0f * I;
        }
    }

    /* previous pilot */

    /*
     * Previous pilot symbol is OFDM_SAMPLESPERFRAME to the left
     * Right after the extra symbol on the left
     *
     * (OFDM_M + OFDM_NCP) + (OFDM_SAMPLESPERFRAME - OFDM_SAMPLESPERFRAME) + 1 + ofdm->sample_point;
     *                                            ^^^^ --- cancels
     */

    st = (OFDM_M + OFDM_NCP) + 1 + ofdm->sample_point;
    en = st + OFDM_M;

    complex float work[OFDM_M];

    /* down-convert at current timing instant---------------------------------- */

    for (j = st, k = 0; j < en; j++, k++) {
        work[k] = ofdm->rxbuf[j] * cexpf(-I * woff_est * j);
    }

    matrix_vector_conjugate_multiply(ofdm, ofdm->rx_sym[0], work); /* sym[0] = previous pilot */

    /* pilot - this frame - pilot */

    /*
     * This pilot is after the extra symbol on the left, and after the OFDM_SAMPLESPERFRAME
     * So we will now be starting at this pilot symbol, and continuing past the next pilot symbol
     * rr = 0 = this pilot, rr = 8 = next pilot
     *
     * Each symbol is of course OFDM_M samples long.
     *
     * In this routine we process the data symbols. These are the ones we want right now.
     */

    for (rr = 0; rr < (OFDM_NS + 1); rr++) {
        st = (OFDM_M + OFDM_NCP) + OFDM_SAMPLESPERFRAME + (rr * (OFDM_M + OFDM_NCP)) + 1 + ofdm->sample_point;
        en = st + OFDM_M;

        /* down-convert at current timing instant---------------------------------- */

        for (j = st, k = 0; j < en; j++, k++) {
            work[k] = ofdm->rxbuf[j] * cexpf(-I * woff_est * j);
        }

        matrix_vector_conjugate_multiply(ofdm, ofdm->rx_sym[rr + 1], work); /* sym[1..9] = this pilot + Data + next pilot */
    }

    /* next pilot */

    /*
     * OK, now we want to process to the next pilot symbol. This is OFDM_SAMPLESPERFRAME symbols to the right
     * Remember now, we are 1 extra symbol + (2 * SAMPLESPERFRAME) symbols from the start right now
     * So in effect, we want to go (3 * SAMPLESPERFRAME) to the right, +/- sample_point
     *
     * We are ignoring the data symbols between the previous pilot and here, we just want the pilot.
     */

    st = (OFDM_M + OFDM_NCP) + (3 * OFDM_SAMPLESPERFRAME) + 1 + ofdm->sample_point;
    en = st + OFDM_M;

    /* down-convert at current timing instant---------------------------------- */

    for (j = st, k = 0; j < en; j++, k++) {
        work[k] = ofdm->rxbuf[j] * cexpf(-I * woff_est * j);
    }

    matrix_vector_conjugate_multiply(ofdm, ofdm->rx_sym[OFDM_NS + 2], work); /* sym[10] = last pilot */

    /*
     * We are finished now with the DFT and down conversion
     * From here on down we are in the frequency domain
     */

    /* est freq err based on all carriers ------------------------------------ */

    if (ofdm->foff_est_en == true) {
        /*
         * Subtract the two complex pilot frequency values
         * sym[1] is the 18 (this) pilot carriers, sym[9] is (next) 18 pilot carriers.
         *
         * So we sum the carriers all together. Since we hope the two pilots are identical,
         * then by conjugating one of them, we should get a real only value, and the
         * imaginary part should be zero.
         *
         * If we find the angle using atan2, then it would be atan2(im, re).
         * but if the im is zero, then we basically at 0 degrees radian at some real x.
         *
         * If the two pilots are different, then we will get a positive or negative angle
         * in radians, and we can use this to correct the frequency.
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

        /*
         * sum the three carrier values, and they should cancel
         * giving us an imaginary value of 0, and thus an angle of 0 degrees radians.
         */

        aphase_est_pilot_rect = aphase_est_pilot_rect + vector_sum(symbol, 3);

        /*
         * Note that i has an index of 1 to 16, so the carriers 0 and 17 will
         * be set to default, or 10 degrees radians, as assigned above.
         *
         * I'm assuming they are just filler, as we don't need them to decode
         * and correct the data symbols. There are only 16 data samples.
         */

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
         * Also note we are using sym[2..8] or the 7 data symbols.
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

