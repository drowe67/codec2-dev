/*
 * Copyright (C) 2017 David Rowe
 *
 * All rights reserved
 *
 * Licensed under GNU LGPL V2.1
 * See LICENSE file for information
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

/*
 * Library of functions that implement a BPSK/QPSK OFDM modem
 *
 * Translated from Octave by Steve Sampson
 */

/* Static Prototypes */

static void matrix_vector_multiply(struct OFDM *, complex float *, complex float *);
static complex float qpsk_mod(int *);
static void qpsk_demod(complex float, int *);
static void ofdm_txframe(struct OFDM *, complex float [OFDM_SAMPLESPERFRAME], complex float *);
static int coarse_sync(struct OFDM *, complex float *, int);

/* Constants */

/*
 * QPSK Quadrant bit-pair values - Gray Coded
 *
 *   0.0 -  89.9 = 00
 *  90.0 - 179.9 = 01
 * 180.0 - 269.9 = 11
 * 270.0 - 359.9 = 10
 */
const complex float constellation[] = {
     1.0f + 0.0f * I,
     0.0f + 1.0f * I,
     0.0f - 1.0f * I,
    -1.0f + 0.0f * I
};

/*
 * These pilots are compatible with Octave version
 */
const char pilotvalues[] = {
    -1, -1, 1, 1, -1, -1, -1, 1, -1,
     1, -1, 1, 1,  1,  1,  1, 1,  1
};

/* Functions */

/* Gray coded QPSK modulation function */

static complex float qpsk_mod(int *bits) {
    return constellation[(bits[0] << 1) + bits[1]];
}

/* Gray coded QPSK demodulation function */

static void qpsk_demod(complex float symbol, int *bits) {
    complex float rotate = symbol * cexpf(I * (M_PI / 4.0f));
    bits[0] = crealf(rotate) < 0.0f;
    bits[1] = cimagf(rotate) < 0.0f;
}

static void matrix_vector_multiply(struct OFDM *ofdm, complex float *result, complex float *vector) {
    int row, col;

    for (row = 0; row < OFDM_M; row++) {
        result[row] = 0.0f + 0.0f * I;

        for (col = 0; col < (OFDM_NC + 2); col++) {
            result[row] += (vector[col] * (ofdm->W[row][col] / (float) OFDM_M)); /* complex result */
        }
    }
}

/*
 * Correlates the OFDM pilot symbol samples with a window of received
 * samples to determine the most likely timing offset.  Combines two
 * frames pilots so we need at least Nsamperframe+M+Ncp samples in rx.
 */

static int coarse_sync(struct OFDM *ofdm, complex float *rx, int length) {
    int Npsam = (OFDM_M + OFDM_NCP);
    int Ncorr = length - (OFDM_SAMPLESPERFRAME + Npsam);
    int i;

    assert(Ncorr > 0);

    float corr[Ncorr];

    for (i = 0; i < Ncorr; i++) {
        complex float csam = conjf(ofdm->rate_fs_pilot_samples[i]);

        corr[i]  = cabsf(rx[i] * csam);
        corr[i] += cabsf(rx[i + OFDM_SAMPLESPERFRAME] * csam);
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
 * ----------------------------------------
 * ofdm_txframe - modulates one frame of symbols
 * ----------------------------------------
 *
 * Each carrier amplitude is 1/M.  There are two edge carriers that
 * are just tx-ed for pilots plus plus Nc continuous carriers. So
 * power is:
 *
 *  p = 2/(Ns*(M*M)) + Nc/(M*M)
 *
 *  e.g. Ns=8, Nc=16, M=144
 *
 * p = 2/(8*(144*144)) + 16/(144*144) = 7.84-04
 *
 */

static void ofdm_txframe(struct OFDM *ofdm, complex float tx[OFDM_SAMPLESPERFRAME],
        complex float *tx_sym_lin) {
    complex float aframe[OFDM_NS][OFDM_NC + 2];
    complex float asymbol[OFDM_M];
    complex float asymbol_cp[OFDM_M + OFDM_NCP];
    int i, j, k, l;

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

    for (i = 0, l = 0; i < OFDM_NS; i++, l += (OFDM_M + OFDM_NCP)) {
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
            tx[l + j] = asymbol_cp[j];
        }
    }
}

/*
 * ------------------------------------------------------------
 * ofdm_create
 *-------------------------------------------------------------
 *
 * Frame has Ns - 1 OFDM data symbols between pilots
 * e.g. for Ns = 3:
 *
 *   PPP
 *   DDD
 *   DDD
 *   PPP
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

    /* store complex pilot symbols */

    for (i = 0; i < (OFDM_NC + 2); i++) {
        ofdm->pilots[i] = ((float) pilotvalues[i]) + 0.0f * I;
    }

    /* carrier tables for up and down conversion */

    int Nlower = floorf((OFDM_CENTRE - OFDM_RS * (OFDM_NC / 2)) / OFDM_RS);

    for (i = 0, j = Nlower; i < (OFDM_NC + 2); i++, j++) {
        ofdm->w[i] = j * TAU * OFDM_RS / OFDM_FS;
    }

    for (i = 0; i < OFDM_M; i++) {
        for (j = 0; j < (OFDM_NC + 2); j++) {
            ofdm->W[i][j] = cexpf(I * ofdm->w[j] * i);
        }
    }

    /* default settings of options and states */

    ofdm->verbose = 0;
    ofdm->timing_en = true;
    ofdm->foff_est_en = true;
    ofdm->phase_est_en = true;

    ofdm->foff_est_gain = 0.01f;
    ofdm->foff_est_hz = 0.0f;
    ofdm->sample_point = 1;
    ofdm->timing_est = 1;
    ofdm->nin = OFDM_SAMPLESPERFRAME;

    /* create the OFDM waveform */

    complex float temp[OFDM_M + OFDM_NCP];

    matrix_vector_multiply(ofdm, temp, ofdm->pilots);

    /*
     * rate_fs_pilot_samples is 160 samples, as we take the last 16 and copy to the front
     * Thus resulting in 16 + 128 + 16 = 160
     */

    /* first copy the last Cyclic Prefix (CP) values */

    for (i = 0, j = (OFDM_M - OFDM_NCP); i < OFDM_NCP; i++, j++) {
        ofdm->rate_fs_pilot_samples[i] = temp[j];
    }

    /* Now copy the whole thing after the above */

    for (i = OFDM_NCP, j = 0; j < OFDM_M; i++, j++) {
        ofdm->rate_fs_pilot_samples[i] = temp[j];
    }

    return ofdm;     /* Success */
}

void ofdm_destroy(struct OFDM *ofdm) {
}

void set_verbose(struct OFDM *ofdm, int level) {
    ofdm->verbose = level;
}

void set_timing_enable(struct OFDM *ofdm, bool val) {
    ofdm->timing_en = val;
}

void set_foff_est_enable(struct OFDM *ofdm, bool val) {
    ofdm->foff_est_en = val;
}

void set_phase_est_enable(struct OFDM *ofdm, bool val) {
    ofdm->phase_est_en = val;
}

void set_foff_est_gain(struct OFDM *ofdm, float val) {
    ofdm->foff_est_gain = val;
}

void set_off_est_hz(struct OFDM *ofdm, float val) {
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
    int s, j;

    if (OFDM_BPS == 1) {
        /* Here we will have Nbitsperframe / 1 */

        for (s = 0; s < length; s++) {
            tx_sym_lin[s] = (float) (2 * tx_bits[s] - 1) + 0.0f * I;
        }
    } else if (OFDM_BPS == 2) {
        /* Here we will have Nbitsperframe / 2 */

        for (s = 0; s < length; s += 2) {
            dibit[0] = tx_bits[s];
            dibit[1] = tx_bits[s + 1];
            tx_sym_lin[s] = qpsk_mod(dibit);
        }
    }

    ofdm_txframe(ofdm, tx, tx_sym_lin);

    /* convert to comp */

    for (s = 0; s < OFDM_SAMPLESPERFRAME; s++) {
        result[s].real = crealf(tx[s]);
        result[s].imag = cimagf(tx[s]);
    }
}

//function [rx_bits states aphase_est_pilot_log rx_np rx_amp] =
int *ofdm_demod(struct OFDM *ofdm, COMP *rxbuf_in) {
/* TODO */
}

