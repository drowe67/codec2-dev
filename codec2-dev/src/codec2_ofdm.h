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

struct OFDM *ofdm_create(float, float, int, float, float, int, int);
void ofdm_destroy(struct OFDM *);
int ofdm_errno(void);
COMP *ofdm_mod(struct OFDM *ofdm, int *);
int *ofdm_demod(struct OFDM *ofdm, COMP *);


#ifdef __cplusplus
}
#endif

#endif
