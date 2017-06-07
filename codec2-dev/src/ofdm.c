/*
 * Copyright (C) 2017 David Rowe
 *
 * All rights reserved
 *
 * Licensed under GNU LGPL V2.1
 * See LICENSE file for information
 */

#define CORRELATE

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <complex.h>

#include "comp.h"
#include "codec2_ofdm.h"
#include "ofdm_internal.h"

#ifdef CORRELATE
#include "codec2_fft.h"
#include "comp_prim.h"

#define FFT_SIZE 512

static codec2_fft_cfg fft_cfg;
#endif

/*
 * Library of functions that implement a BPSK/QPSK OFDM modem
 *
 * Translated from Octave by Steve Sampson
 */

/* Static Prototypes */

static void matrix_vector_multiply(complex float *, int, int, complex float **, int *);
static complex float qpsk_mod(int *);
static void qpsk_demod(complex float, int *);
static void ofdm_txframe(struct OFDM *, complex float **, complex float *);
static int coarse_sync(struct OFDM *, int *, complex float *, int);

/* global error value */

static int errno;

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

static void matrix_vector_multiply(complex float *result,
        int columns, int rows, complex float **array, int *vector) {
    int i, j, k;

    for (i = 0; i < columns; i++) {
        for (j = 0; j < rows; j++) {
            result[j] = CMPLXF(0.0f, 0.0f);

            for (k = 0; k < columns; k++) {
                result[j] += (vector[k] * (array[k][j] / (float) rows));    /* complex result */
            }
        }
    }
}

/*
 * Correlates the OFDM pilot symbol samples with a window of received
 * samples to determine the most likely timing offset.  Combines two
 * frames pilots so we need at least Nsamperframe+M+Ncp samples in rx.
 *
 * Also determines frequency offset at maximum correlation.  Can be
 * used for acquisition (coarse timing a freq offset), and fine timing
 */

static int coarse_sync(struct OFDM *ofdm, int *foff_est, complex float *rx, int length) {
    int Npsam = (ofdm->M + ofdm->Ncp);
    int Ncorr = length - (ofdm->Nsamperframe + Npsam);
    int i;

    assert(Ncorr > 0);

    float corr[Ncorr];

    for (i = 0; i < Ncorr; i++) {
        complex float csam = conjf(ofdm->rate_fs_pilot_samples[i]);

        corr[i]  = cabsf(rx[i] * csam);
        corr[i] += cabsf(rx[i + ofdm->Nsamperframe] * csam);
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

#ifdef CORRELATE
    COMP binl[FFT_SIZE];
    COMP binr[FFT_SIZE];
    int i, j;

    /* determines frequency offset at maximum correlation */

    for (i = t_est, j = 0; i < (t_est + Npsam); i++, j++) {
        complex float csam = conjf(ofdm->rate_fs_pilot_samples[j]);
        complex float templ = (rx[i] * csam);
        complex float tempr = (rx[i + ofdm->Nsamperframe] * csam);

        binl[j].real = crealf(templ);
        binl[j].imag = cimagf(templ);

        binr[j].real = crealf(tempr);
        binr[j].imag = cimagf(tempr);
    }

    /* Initialize rest of the FFT bins */

    for (i = Npsam; i < FFT_SIZE; i++) {
        binl[i].real = 0.0f;
        binl[i].imag = 0.0f;
        
        binr[i].real = 0.0f;
        binr[i].imag = 0.0f;
    }

    codec2_fft_inplace(fft_cfg, binl);
    codec2_fft_inplace(fft_cfg, binr);

    int fmax = 30;
    float C[fmax] = 0.0f;

    for (i = 0; i < (FFT_SIZE / 2); i++) {
        C[i] = cabsolute(binl[i]);
        C[i] += cabsolute(binr[i]);
    }

    float mx_pos = 0.0f;
    float mx_neg = 0.0f;

    int foff_est_pos = 0;
    int foff_est_neg = 0;

    for (i = 0; i < fmax; i++) {
        if (C[i] > mx_pos) {
            mx_pos = C[i];
            foff_est_pos = i;
        }
    }

    for (i = ((FFT_SIZE / 2) - fmax); i < (FFT_SIZE / 2); i++) {
        if (C[i] > mx_neg) {
            mx_neg = C[i];
            foff_est_neg = i;
        }
    }

    if (mx_pos > mx_neg) {
      *foff_est = foff_est_pos;
    } else {
      *foff_est = foff_est_neg - fmax;
    }
#endif

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

static void ofdm_txframe(struct OFDM *ofdm, complex float **tx, complex float *tx_sym_lin) {
    complex float aframe[ofdm->Nrowsperframe][ofdm->Nc + 2];   /* [8][18] */
    int i, j, k;

    /* initialize aframe to complex zero */

    for (i = 0; i < ofdm->Nrowsperframe; i++) {
        for (j = 0; j < (ofdm->Nc + 2); j++) {
            aframe[i][j] = CMPLXF(0.0f, 0.0f);
        }
    }

    /* copy in a row of complex pilots to first row */

    for (i = 0; i < (ofdm->Nc + 2); i++) {
        aframe[0][i] = ofdm->pilots[i];
    }

    /* Place symbols in multi-carrier frame with pilots */
    /* This will place boundary values of complex zero around data */

    for (i = 1; i < ofdm->Nrowsperframe; i++) {

        /* copy in the Nc complex values with [0 Nc 0] or (Nc + 2) total */


        for (j = 1; j < (ofdm->Nc + 1); j++) {
            aframe[i][j] = tx_sym_lin[((i - 1) * ofdm->Nc) + (j - 1)];
        }
    }

    /* OFDM up-convert symbol by symbol so we can add CP */

    for (i = 0; i < ofdm->Nrowsperframe; i++) {
        complex float asymbol[ofdm->M];
        complex float asymbol_cp[ofdm->M + ofdm->Ncp];

        matrix_vector_multiply(&asymbol, (ofdm->Nc + 2), ofdm->M, ofdm->W, &aframe[i]);

        /* Copy the last Ncp + 2 columns to the front */

        for (j = (ofdm->M - (ofdm->Ncp + 2)), k = 0; j < ofdm->M; j++, k++) {
            asymbol_cp[k] = asymbol[j];
        }

        /* Now copy the whole row after it */

        for (j = (ofdm->Ncp + 2), k = 0; j < ofdm->M; j++, k++) {
            asymbol_cp[j] = asymbol[k];
        }

        /* Now move row to the tx reference */

        for (j = 0; j < (ofdm->M + ofdm->Ncp); j++) {
            tx[i][j] = asymbol_cp[j];
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
 * Return NULL with errno set to error code on fail
 */

struct OFDM *ofdm_create(float freq, float fs, int bps, float ts, float tcp, int ns, int nc) {
    struct OFDM *ofdm;
    int i, j;

    if ((ofdm = (struct OFDM *) malloc(sizeof (struct OFDM))) == NULL) {
	return NULL;
    }

    ofdm->Fcentre = freq;
    ofdm->Fs = fs;                        /* Sample Rate */
    ofdm->bps = bps;                      /* bits per symbol (1 = BPSK/ 2 = QPSK) */
    ofdm->Ts = ts;                        /* Symbol Time .018 */
    ofdm->Rs = (1.0f / ofdm->Ts);         /* Symbol Rate */
    ofdm->Tcp = tcp;
    ofdm->Ns = ns;
    ofdm->Nc = nc;
    ofdm->M = (int)(ofdm->Fs / ofdm->Rs);                                   /* 144 */
    ofdm->Ncp = (int)(ofdm->Tcp * ofdm->Fs);                                /* 16 */
    ofdm->Nbitsperframe = (ofdm->Ns * ofdm->Nc * ofdm->bps);                /* 256 */
    ofdm->Nrowsperframe = (ofdm->Nbitsperframe / (ofdm->Nc * ofdm->bps));   /* 8 */
    ofdm->Nsamperframe = (ofdm->Nrowsperframe * (ofdm->M + ofdm->Ncp));     /* 1280 */

    if (ofdm->Fs < 8000.0f) {
        errno = 1;
        return NULL;
    }

    if ((ofdm->bps < 1) || (ofdm->bps > 2)) {
        errno = 2;
        return NULL;
    }

    if (ofdm->Ts >= 1.0f) {
        errno = 3;
        return NULL;
    }

    /*
     * Since the modem can be configured using different parameters, we use dynamic
     * memory allocation to build the data structures.
     *
     * Whether malloc or calloc is called, is dependent on whether the array or vector
     * is assigned values here that completely fill the structure.
     */

    if ((ofdm->rate_fs_pilot_samples = (complex float *) malloc((ofdm->M + ofdm->Nc) * sizeof (complex float))) == NULL) {
        errno = 1000;
        return NULL;
    }

    if ((ofdm->w = (float *) malloc((ofdm->Nc + 2) * sizeof (float))) == NULL) {
        free(ofdm->rate_fs_pilot_samples);
        errno = 2000;
        return NULL;
    }

    /*
     * Receive buffer: D P DDD P DDD P DDD P D
     *                         ^
     * also see ofdm_demod() ...
     *
     * 4320 complex floats (@ 34560 bytes)
     *
     * This buffer should be zero'd out as per octave
     */

    ofdm->Nrxbuf = 3 * ofdm->Nsamperframe + 3 * (ofdm->M + ofdm->Ncp);

    if ((ofdm->rxbuf = (complex float *) calloc(ofdm->Nrxbuf, sizeof (complex float))) == NULL) {
        free(ofdm->w);
        free(ofdm->rate_fs_pilot_samples);
        errno = 3000;
        return NULL;
    }

    if ((ofdm->pilots = (int *) malloc((ofdm->Nc + 2) * sizeof (int))) == NULL) {
        free(ofdm->rxbuf);
        free(ofdm->w);
        free(ofdm->rate_fs_pilot_samples);
        errno = 4000;
        return NULL;
    }

    /* Allocate Nc + 2 columns of complex pointers */

    if ((ofdm->W = (complex float **) malloc((ofdm->Nc + 2) * sizeof (complex float *))) == NULL) {
        free(ofdm->pilots);
        free(ofdm->rxbuf);
        free(ofdm->w);
        free(ofdm->rate_fs_pilot_samples);
        errno = 5000;
        return NULL;
    }

    /* Allocate M rows of complex pointers */

    for (i = 0; i < ofdm->M; i++) {
        if ((*(ofdm->W + i) = (complex float *) malloc(ofdm->M * sizeof (complex float))) == NULL) {
            free(ofdm->W);
            free(ofdm->pilots);
            free(ofdm->rxbuf);
            free(ofdm->w);
            free(ofdm->rate_fs_pilot_samples);
            errno = 6000;
            return NULL;
        }
    }

#ifdef CORRELATE
    if ((fft_cfg = codec2_fft_alloc(FFT_SIZE, 0, NULL, NULL)) == NULL) {
        for (i = 0; i < ofdm->M; i++) {
            free(*(ofdm->W + i));
        }

        free(ofdm->W);
        free(ofdm->pilots);
        free(ofdm->rxbuf);
        free(ofdm->w);
        free(ofdm->rate_fs_pilot_samples);
        errno = 7000;
        return NULL;
    }
#endif

    /* same pilots each time */

    srand(1);

    /* store pilot symbols in allocated memory */

    for (i = 0; i < (ofdm->Nc + 2); i++) {
        ofdm->pilots[i] = 1 - 2 * (((float)rand() / (float)RAND_MAX) > 0.5f);
    }

    /* carrier tables for up and down conversion */

    int Nlower = floorf((ofdm->Fcentre - ofdm->Rs * (ofdm->Nc / 2)) / ofdm->Rs);

    for (i = 0, j = Nlower; i < (ofdm->Nc + 2); i++, j++) {
        ofdm->w[i] = j * TAU * ofdm->Rs / ofdm->Fs;
    }

    for (i = 0; i < (ofdm->Nc + 2); i++) {
        for (j = 0; j < ofdm->M; j++) {
            ofdm->W[i][j] = cexpf(I * ofdm->w[i] * j);
        }
    }

    /* fine timing search +/- window_width/2 from current timing instant */

    ofdm->ftwindow_width = 11;

    /* default settings of options and states */

    ofdm->verbose = 0;
    ofdm->timing_en = true;
    ofdm->foff_est_en = true;
    ofdm->phase_est_en = true;

    ofdm->foff_est_gain = 0.01f;
    ofdm->foff_est_hz = 0.0f;
    ofdm->sample_point = 1;
    ofdm->timing_est = 1;
    ofdm->nin = ofdm->Nsamperframe;

    /* create the OFDM waveform */

    complex float temp[ofdm->Ncp + ofdm->M];

    matrix_vector_multiply(temp, (ofdm->Nc + 2), ofdm->M, ofdm->W, ofdm->pilots);

    /*
     * rate_fs_pilot_samples is 160 samples, as we take the last 16 and copy to the front
     * Thus resulting in 16 + 128 + 16 = 160
     */

    /* first copy the last Ncp values */

    for (i = 0, j = (ofdm->M - ofdm->Ncp); i < (ofdm->M - ofdm->Ncp); i++, j++) {
        ofdm->rate_fs_pilot_samples[i] = temp[j];
    }

    /* Now copy the whole thing after the above */

    for (i = ofdm->Ncp, j = 0; i < ofdm->M; i++, j++) {
        ofdm->rate_fs_pilot_samples[i] = temp[j];
    }

    errno = 0;

    return ofdm;     /* Success */
}

void ofdm_destroy(struct OFDM *ofdm) {
    int i;

#ifdef CORRELATE
    codec2_fft_free(fft_cfg);
#endif
    for (i = 0; i < ofdm->M; i++) {
        free(*(ofdm->W + i));
    }

    free(ofdm->W);
    free(ofdm->pilots);
    free(ofdm->rxbuf);
    free(ofdm->w);
    free(ofdm->rate_fs_pilot_samples);
}

/*
 * Optional function in case the user wants
 * to know where the library failed.
 */

int ofdm_errno() {
    return errno;
}

void set_verbose(struct OFDM *ofdm, int level) {
    ofdm->verbose = level;
}

int get_verbose(struct OFDM *ofdm) {
    return ofdm->verbose;
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

COMP *ofdm_mod(struct OFDM *ofdm, int *tx_bits) {
    int length = ofdm->Nbitsperframe / ofdm->bps;
    complex float tx[ofdm->Nrowsperframe][ofdm->M + ofdm->Ncp];
    complex float tx_sym_lin[length];
    struct COMP result[ofdm->Nrowsperframe][ofdm->M + ofdm->Ncp];
    int i, j;

    switch (ofdm->bps) {
        case 1:
            /* Here we will have Nbitsperframe / 1 */

            for (int s = 0; s < length; s++) {
                tx_sym_lin[s] = CMPLXF((float)(2 * tx_bits[s] - 1), 0.0f);
            }
            break;
        default:
        case 2:
            /* Here we will have Nbitsperframe / 2 */

            for (int s = 0; s < length; s++) {
                tx_sym_lin[s] = qpsk_mod(tx_bits);
            }
    }

    ofdm_txframe(ofdm, tx, tx_sym_lin);

    /* convert to comp */

    for (i = 0; i < ofdm->Nrowsperframe; i++) {
        for (j = 0; j < (ofdm->M + ofdm->Ncp); j++) {
            result[i][j].real = crealf(tx[i][j]);
            result[i][j].imag = cimagf(tx[i][j]);
        }
    }

    return result;
}

/*
 * ------------------------------------------
 * ofdm_demod - Demodulates one frame of bits
 * ------------------------------------------
 *
 * For phase estimation we need to maintain buffer of 3 frames plus
 * one pilot, so we have 4 pilots total. '^' is the start of current
 * frame that we are demodulating.
 *
 * P DDD P DDD P DDD P
 *       ^
 *
 * Then add one symbol either side to account for movement in
 * sampling instant due to sample clock differences:
 *
 * D P DDD P DDD P DDD P D
 *         ^
 */

//function [rx_bits states aphase_est_pilot_log rx_np rx_amp] =
int *ofdm_demod(struct OFDM *ofdm, COMP *rxbuf_in) {
/* TODO */
}

