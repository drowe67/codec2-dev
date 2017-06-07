/*
 * Copyright (C) 2017 David Rowe
 *
 * All rights reserved
 * 
 * Licensed under GNU LGPL V2.1
 * See LICENSE file for information
 */

#ifndef OFDM_INTERNAL_H
#define OFDM_INTERNAL_H

#ifdef __cplusplus
extern "C" {
#endif

#include <complex.h>
#include <stdbool.h>

#ifndef M_PI
#define M_PI       3.14159265358979323846f  /* math constant */
#endif

#define TAU        (2.0f * M_PI)            /* mathematical constant */

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

struct OFDM {
    float Fs;
    float Ts;
    float Rs;
    float Tcp;
    float Fcentre;
    float foff_est_gain;
    float foff_est_hz;

    int Nbitsperframe;
    int Nrowsperframe;
    int Nsamperframe;
    int Ns;
    int Nc;
    int M;
    int bps;
    int ftwindow_width;
    int verbose;
    int sample_point;
    int timing_est;
    int nin;
    int Nrxbuf;
    int Ncp;

    bool timing_en;
    bool foff_est_en;
    bool phase_est_en;

    /* dynamic heap memory allocation */

    complex float *rate_fs_pilot_samples;
    complex float **W;
    complex float *rxbuf;
    complex float *pilots;
    float *w;
};

#ifdef __cplusplus
}
#endif

#endif
