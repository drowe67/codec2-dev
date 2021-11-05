/*---------------------------------------------------------------------------*\

  FILE........: fsk.h
  AUTHOR......: Brady O'Brien
  DATE CREATED: 6 January 2016

  C Implementation of 2FSK/4FSK modulator/demodulator, based on octave/fsk_horus.m

\*---------------------------------------------------------------------------*/

/*
  Copyright (C) 2016 David Rowe

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


#ifndef __C2FSK_H
#define __C2FSK_H

#define FFTW3           // comment out to use KISS FFT

#include <stdint.h>
#include "comp.h"

#ifdef FFTW3
#include <complex.h>
#include <fftw3.h>
#else
#include "kiss_fftr.h"
#endif

#include "modem_stats.h"

#define MODE_2FSK 2
#define MODE_4FSK 4

#define MODE_M_MAX 4

#define FSK_SCALE 16383

/* default internal parameters */
#define FSK_DEFAULT_P     8      /* Number of timing offsets we have to choose from, try to keep P >= 8 */
#define FSK_DEFAULT_NSYM 50      /* See Nsym below */
#define FSK_NONE         -1      /* unused parameter */

struct FSK {
#ifdef FFTW3
    fftwf_plan plan;
    fftwf_complex *fftin;
    fftwf_complex *fftout;
#endif
    
    /*  Static parameters set up by fsk_init */
    int Ndft;               /* freq offset est fft */
    int Fs;                 /* sample freq */
    int N;                  /* processing buffer size */
    int Rs;                 /* symbol rate */
    int Ts;                 /* samples per symbol */
    int Nmem;               /* size of extra mem for timing adj */
    int P;                  /* oversample rate for timing est/adj */
    int Nsym;               /* Number of symbols processed by demodulator in each call, also the timing estimator window */
    int Nbits;              /* Number of bits spat out in a processing frame */
    int f1_tx;              /* f1 for modulator */
    int tone_spacing;       /* Space between TX freqs for modulator (and option mask freq estimator) */
    int mode;               /* 2FSK or 4FSK */
    float tc;               /* time constant for smoothing FFTs */
    int est_min;            /* Minimum frequency for freq. estimator */
    int est_max;            /* Maximum frequency for freq. estimator */
    int est_space;          /* Minimum frequency spacing for freq. estimator */
    float* hann_table;	    /* Precomputed or runtime computed hann window table */
    
    /*  Parameters used by demod */
    float* Sf;	            /* Average of magnitude spectrum */
    COMP phi_c[MODE_M_MAX]; /* phase of each demod local oscillator */
    COMP *f_dc;             /* down converted samples               */
        
#ifndef FFTW3
    kiss_fft_cfg fft_cfg;   /* Config for KISS FFT, used in freq est */
#endif

    float norm_rx_timing;   /* Normalized RX timing */
            
    /*  Parameters used by mod */
    COMP tx_phase_c;        /* TX phase, but complex */ 
    
    /*  Statistics generated by demod */
    float  EbNodB;            /* Estimated EbNo in dB */
    float  f_est[MODE_M_MAX]; /* Estimated frequencies (peak method) */
    float  f2_est[MODE_M_MAX];/* Estimated frequencies (mask method) */
    int    freq_est_type;     /* which estimator to use              */
    float  ppm;               /* Estimated PPM clock offset */
    float  SNRest;            /* used for LLRs */
    float  v_est;             /* used for LLRs */
    float  rx_sig_pow;
    float  rx_nse_pow;
    
    /*  Parameters used by mod/demod and driving code */
    int nin;                /* Number of samples to feed the next demod cycle */
    int burst_mode;         /* enables/disables 'burst' mode */
    int lock_nin;           /* locks nin during testing */
    
    /*  modem statistics struct */
    struct MODEM_STATS *stats;
    int normalise_eye;      /* enables/disables normalisation of eye diagram */
};

/*
 * Create a FSK modem
 * 
 * int Fs - Sample frequency
 * int Rs - Symbol rate
 * int M  - 2 for 2FSK, 4 for 4FSK
 * int f1_tx - first tone frequency
 * int tone_spacing - frequency spacing (for modulator and optional "mask" freq estimator)
 */
struct FSK * fsk_create(int Fs, int Rs, int M, int f1_tx, int tone_spacing);

/*
 * Create a FSK modem - advanced version
 * 
 * int Fs - Sample frequency
 * int Rs - Symbol rate
 * int M  - 2 for 2FSK, 4 for 4FSK
 * int P  - number of timing offsets to choose from (suggest >= 8)
 * int Nsym  - windows size for timing estimator
 * int f1_tx - first tone frequency
 * int tone_spacing - frequency spacing (for modulator and optional "mask" freq estimator)
 */
struct FSK * fsk_create_hbr(int Fs, int Rs, int M, int P, int Nsym, int f1_tx, int tone_spacing);

/*
 * Set the minimum and maximum frequencies at which the freq. estimator can find tones
 */
void fsk_set_freq_est_limits(struct FSK *fsk,int fmin, int fmax);

/* 
 * Clear the estimator states
 */
void fsk_clear_estimators(struct FSK *fsk);

/*
 * Fills MODEM_STATS struct with demod statistics
 */
void fsk_get_demod_stats(struct FSK *fsk, struct MODEM_STATS *stats);

/*
 * Destroy an FSK state struct and free it's memory
 * 
 * struct FSK *fsk - FSK config/state struct to be destroyed
 */
void fsk_destroy(struct FSK *fsk);

/*
 * Modulates Nsym bits into N samples
 * 
 * struct FSK *fsk   - FSK config/state struct, set up by fsk_create
 * float fsk_out[]   - Buffer for samples of modulated FSK, fsk->Ts*(nbits/(M>>1)) in length
 * uint8_t tx_bits[] - Buffer containing Nbits unpacked bits
 * int     nbits     - number of bits to transmit
 */
void fsk_mod(struct FSK *fsk, float fsk_out[], uint8_t tx_bits[], int nbits);

/*
 * Modulates Nsym bits into N samples
 * 
 * struct FSK *fsk   - FSK config/state struct, set up by fsk_create
 * float fsk_out[]   - Buffer for samples of "voltage" used to modulate an external VCO
 *                   - fsk->Ts*(nbits/(M>>1)) in length
 * uint8_t tx_bits[] - Buffer containing Nbits unpacked bits
 * int     nbits     - number of bits to transmit
 */
void fsk_mod_ext_vco(struct FSK *fsk, float vco_out[], uint8_t tx_bits[], int nbits);

/*
 * Modulates Nsym bits into N complex samples
 * 
 * struct FSK *fsk   - FSK config/state struct, set up by fsk_create
 * comp fsk_out[]    - Buffer for samples of modulated FSK, fsk->Ts*(nbits/(M>>1)) in length
 * uint8_t tx_bits[] - Buffer containing Nbits unpacked bits
 * int     nbits     - number of bits to transmit
 */
void fsk_mod_c(struct FSK *fsk, COMP fsk_out[], uint8_t tx_bits[], int nbits);

/*
 * Returns the number of samples needed for the next fsk_demod() cycle
 *
 * struct FSK *fsk - FSK config/state struct, set up by fsk_create
 * returns - number of samples to be fed into fsk_demod next cycle 
 */
uint32_t fsk_nin(struct FSK *fsk);


/*
 * Demodulate some number of FSK samples. The number of samples to be 
 *  demodulated can be found by calling fsk_nin().
 * 
 * struct FSK *fsk   - FSK config/state struct, set up by fsk_create
 * uint8_t rx_bits[] - Buffer for fsk->Nbits unpacked bits to be written
 * float fsk_in[]    - nin samples of modulated FSK
 */
void fsk_demod(struct FSK *fsk, uint8_t rx_bits[], COMP fsk_in[]);

/*
 * Soft decision demodulation
 * 
 * struct FSK *fsk - FSK config/state struct, set up by fsk_create
 * float rx_flit[] - M x Nsym array of filtermagnitude outputs
 * float fsk_in[]  - nin samples of modualted FSK
 */
void fsk_demod_sd(struct FSK *fsk, float rx_filt[], COMP fsk_in[]);

/* enables/disables normalisation of eye diagram samples */
  
void fsk_stats_normalise_eye(struct FSK *fsk, int normalise_enable);

/* Set the FSK modem into burst demod mode */

void fsk_enable_burst_mode(struct FSK *fsk);

/* Set freq est algorithm 0: peak 1:mask */
void fsk_set_freq_est_alg(struct FSK *fsk, int est_type);

#endif
