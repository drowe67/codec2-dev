/*
 * Copyright (C) 2017 David Rowe
 *
 * All rights reserved
 * 
 * Licensed under GNU LGPL V2.1
 * See LICENSE file for information
 */

#ifndef CODEC2_OFDM_H
#define CODEC2_OFDM_H

#ifdef __cplusplus
extern "C" {
#endif

/* Includes */
    
#include <complex.h>
#include <stdbool.h>
    
#include "comp.h"

/* Defines */

struct OFDM;

/* Prototypes */

struct OFDM *ofdm_create(void);
void ofdm_destroy(struct OFDM *);
void ofdm_mod(struct OFDM *ofdm, COMP [OFDM_SAMPLESPERFRAME], int *);
int *ofdm_demod(struct OFDM *ofdm, COMP *);

/* option setters */

void set_verbose(struct OFDM *, int);
void set_timing_enable(struct OFDM *, bool);
void set_foff_est_enable(struct OFDM *, bool);
void set_phase_est_enable(struct OFDM *, bool);
void set_foff_est_gain(struct OFDM *, float);
void set_off_est_hz(struct OFDM *, float);

#ifdef __cplusplus
}
#endif

#endif

