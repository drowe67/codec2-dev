/*---------------------------------------------------------------------------*\

  FILE........: ofdm.c
  AUTHORS.....: David Rowe & Steve Sampson
  DATE CREATED: June 2017

  A Library of functions that implement a QPSK OFDM modem, C port of
  the Octave functions in ofdm_lib.m

\*---------------------------------------------------------------------------*/
/*
  Copyright (C) 2017-2019 David Rowe

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
#include "filter.h"
#include "wval.h"
#include "debug_alloc.h"
#include "machdep.h"

/* Static Prototypes */

static float cnormf(complex float);
static void allocate_tx_bpf(struct OFDM *);
static void deallocate_tx_bpf(struct OFDM *);
static void dft(struct OFDM *, complex float *, complex float *);
static void idft(struct OFDM *, complex float *, complex float *);
static complex float vector_sum(complex float *, int);
static int est_timing(struct OFDM *, complex float *, int, int, float *, int *, int);
static float est_freq_offset_pilot_corr(struct OFDM *, complex float *, int, int);
static int ofdm_sync_search_core(struct OFDM *);
static void ofdm_demod_core(struct OFDM *, int *);

/* Defines */

#define max( a, b ) ( ((a) > (b)) ? (a) : (b) )
#define min( a, b ) ( ((a) < (b)) ? (a) : (b) )

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
static const int8_t pilotvalues[] = {
  -1,-1, 1, 1,-1,-1,-1, 1,
  -1, 1,-1, 1, 1, 1, 1, 1,
   1, 1, 1,-1,-1, 1,-1, 1,
  -1, 1, 1, 1, 1, 1, 1, 1,
   1, 1, 1,-1, 1, 1, 1, 1,
   1,-1,-1,-1,-1,-1,-1, 1,
  -1, 1,-1, 1,-1,-1, 1,-1,
   1, 1, 1, 1,-1, 1,-1, 1
};

/* static variables - used as storage for constants that are unchanged after init-time */

static struct OFDM_CONFIG ofdm_config;

static complex float *tx_uw_syms;
static int *uw_ind;
static int *uw_ind_sym;

static float ofdm_tx_centre; /* TX Center frequency */
static float ofdm_rx_centre; /* RX Center frequency */
static float ofdm_fs; /* Sample rate */
static float ofdm_ts; /* Symbol cycle time */
static float ofdm_rs; /* Symbol rate */
static float ofdm_tcp; /* Cyclic prefix duration */
static float ofdm_inv_m; /* 1/m */
static float ofdm_tx_nlower; /* TX lowest carrier freq */
static float ofdm_rx_nlower; /* RX lowest carrier freq */
static float ofdm_doc; /* division of radian circle */

/*
 * See 700D Part 4 Acquisition blog post and ofdm_dev.m routines
 * for how this was set
 */
static float ofdm_timing_mx_thresh;

static int ofdm_nc;
static int ofdm_ns;	/* NS-1 = data symbols between pilots  */
static int ofdm_bps; 	/* Bits per symbol */
static int ofdm_m; 	/* duration of each symbol in samples */
static int ofdm_ncp; 	/* duration of CP in samples */

static int ofdm_ftwindowwidth;
static int ofdm_bitsperframe;
static int ofdm_rowsperframe;
static int ofdm_samplesperframe;
static int ofdm_max_samplesperframe;
static int ofdm_rxbuf;
static int ofdm_ntxtbits; /* reserve bits/frame for aux text information */
static int ofdm_nuwbits; /* Unique word used for positive indication of lock */

/* Local Functions ----------------------------------------------------------*/

static float cnormf(complex float val) {
    float realf = crealf(val);
    float imagf = cimagf(val);

    return realf * realf + imagf * imagf;
}

/*
 * Gray coded QPSK modulation function
 */
complex float qpsk_mod(int *bits) {
    return constellation[(bits[1] << 1) | bits[0]];
}

/*
 * Gray coded QPSK demodulation function
 *
 * 01 | 00
 * ---+---
 * 11 | 10
 */
void qpsk_demod(complex float symbol, int *bits) {
    complex float rotate = symbol * cmplx(ROT45);

    bits[0] = crealf(rotate) < 0.0f;
    bits[1] = cimagf(rotate) < 0.0f;
}

/*
 * ------------
 * ofdm_create
 * ------------
 *
 * Returns OFDM data structure on success
 * Return NULL on fail
 *
 * If you want the defaults, call this with config structure
 * and the NC setting to 0. This will fill the structure with
 * default values of the original OFDM modem.
 */
struct OFDM *ofdm_create(const struct OFDM_CONFIG *config) {
    struct OFDM *ofdm;
    float tval;
    int i, j;

    /* Check if called correctly */

    if (config == NULL) {
        return NULL;
    }

    if (config->nc == 0) {
        /* Fill in default values */

        ofdm_nc = 17; /* Number of carriers */
        ofdm_ns = 8; /* Number of Symbol frames */
        ofdm_bps = 2; /* Bits per Symbol */
        ofdm_ts = 0.018f;
        ofdm_tcp = .002f; /* Cyclic Prefix duration */
        ofdm_tx_centre = 1500.0f; /* TX Centre Audio Frequency */
        ofdm_rx_centre = 1500.0f; /* RX Centre Audio Frequency */
        ofdm_fs = 8000.0f; /* Sample Frequency */
        ofdm_ntxtbits = 4;
        ofdm_ftwindowwidth = 11;
        ofdm_timing_mx_thresh = 0.30f;
    } else {
        /* Use the users values */

        ofdm_nc = config->nc; /* Number of carriers */
        ofdm_ns = config->ns; /* Number of Symbol frames */

        if ((config->bps != 1) && (config->bps != 2)) {
            ofdm_bps = 2;   /* punt on bad data */
        } else {
            ofdm_bps = config->bps; /* Bits per Symbol */
        }

        ofdm_ts = config->ts;
        ofdm_tcp = config->tcp; /* Cyclic Prefix duration */
        ofdm_tx_centre = config->tx_centre; /* TX Centre Audio Frequency */
        ofdm_rx_centre = config->rx_centre; /* RX Centre Audio Frequency */
        ofdm_fs = config->fs; /* Sample Frequency */
        ofdm_ntxtbits = config->txtbits;

        ofdm_ftwindowwidth = config->ftwindowwidth;
        ofdm_timing_mx_thresh = config->ofdm_timing_mx_thresh;
    }

    ofdm_rs = (1.0f / ofdm_ts); /* Modulation Symbol Rate */
    ofdm_m = (int) (ofdm_fs / ofdm_rs); /* 144 */
    ofdm_ncp = (int) (ofdm_tcp * ofdm_fs); /* 16 */
    ofdm_inv_m = (1.0f / (float) ofdm_m);

    /* Copy structure into global */

    ofdm_config.tx_centre = ofdm_tx_centre;
    ofdm_config.rx_centre = ofdm_rx_centre;
    ofdm_config.fs = ofdm_fs;
    ofdm_config.rs = ofdm_rs;
    ofdm_config.ts = ofdm_ts;
    ofdm_config.tcp = ofdm_tcp;
    ofdm_config.ofdm_timing_mx_thresh = ofdm_timing_mx_thresh;
    ofdm_config.nc = ofdm_nc;
    ofdm_config.ns = ofdm_ns;
    ofdm_config.bps = ofdm_bps;
    ofdm_config.txtbits = ofdm_ntxtbits;
    ofdm_config.ftwindowwidth = ofdm_ftwindowwidth;

    /* Calculate sizes from config param */

    ofdm_bitsperframe = (ofdm_ns - 1) * (ofdm_nc * ofdm_bps);
    ofdm_rowsperframe = ofdm_bitsperframe / (ofdm_nc * ofdm_bps);
    ofdm_samplesperframe = ofdm_ns * (ofdm_m + ofdm_ncp);
    ofdm_max_samplesperframe = ofdm_samplesperframe + (ofdm_m + ofdm_ncp) / 4;
    ofdm_rxbuf = 3 * ofdm_samplesperframe + 3 * (ofdm_m + ofdm_ncp);
    ofdm_nuwbits = (ofdm_ns - 1) * ofdm_bps - ofdm_ntxtbits;    // 10
    
    /* Were ready to start filling in the OFDM structure now */
    ofdm = (struct OFDM *) MALLOC(sizeof (struct OFDM));
    assert(ofdm != NULL);

    ofdm->pilot_samples = MALLOC(sizeof (complex float) * (ofdm_m + ofdm_ncp));
    assert(ofdm->pilot_samples != NULL);

    ofdm->rxbuf = MALLOC(sizeof (complex float) * ofdm_rxbuf);
    assert(ofdm->rxbuf != NULL);

    ofdm->pilots = MALLOC(sizeof (complex float) * (ofdm_nc + 2));
    assert(ofdm->pilots !=  NULL);

    /*
     * rx_sym is a 2D array of variable size
     *
     * allocate rx_sym row storage. It is a pointer to a pointer
     */
    ofdm->rx_sym = MALLOC(sizeof (complex float) * (ofdm_ns + 3));
    assert(ofdm->rx_sym != NULL);

    /* allocate rx_sym column storage */

    for (i = 0; i < (ofdm_ns + 3); i++) {
        ofdm->rx_sym[i] = (complex float *) MALLOC(sizeof(complex float) * (ofdm_nc + 2));
	assert(ofdm->rx_sym[i] != NULL);
    }

    /* The rest of these are 1D arrays of variable size */

    ofdm->rx_np = MALLOC(sizeof (complex float) * (ofdm_rowsperframe * ofdm_nc));
    assert(ofdm->rx_np != NULL);

    ofdm->rx_amp = MALLOC(sizeof (float) * (ofdm_rowsperframe * ofdm_nc));
    assert(ofdm->rx_amp != NULL);

    ofdm->aphase_est_pilot_log = MALLOC(sizeof (float) * (ofdm_rowsperframe * ofdm_nc));
    assert(ofdm->aphase_est_pilot_log != NULL);

    ofdm->tx_uw = MALLOC(sizeof (uint8_t) * ofdm_nuwbits);
    assert(ofdm->tx_uw != NULL);

    for (i = 0; i < ofdm_nuwbits; i++) {
        ofdm->tx_uw[i] = 0;
    }

    /* Null pointers to unallocated buffers */
    ofdm->ofdm_tx_bpf = NULL;

    /* store complex BPSK pilot symbols */

    assert(sizeof (pilotvalues) >= (ofdm_nc + 2) * sizeof (int8_t));

    /* There are only 64 pilot values available */

    for (i = 0; i < (ofdm_nc + 2); i++) {
        ofdm->pilots[i] = ((float) pilotvalues[i]) + 0.0f * I;
    }

    /* carrier tables for up and down conversion */

    ofdm_doc = (TAU / (ofdm_fs / ofdm_rs));
    tval = ((float) ofdm_nc / 2);
    ofdm_tx_nlower = roundf((ofdm_tx_centre / ofdm_rs) - tval) - 1;
    ofdm_rx_nlower = roundf((ofdm_rx_centre / ofdm_rs) - tval) - 1;

    for (i = 0; i < ofdm_rxbuf; i++) {
        ofdm->rxbuf[i] = 0.0f;
    }

    for (i = 0; i < (ofdm_ns + 3); i++) {
        for (j = 0; j < (ofdm_nc + 2); j++) {
            ofdm->rx_sym[i][j] = 0.0f;
        }
    }

    for (i = 0; i < ofdm_rowsperframe * ofdm_nc; i++) {
        ofdm->rx_np[i] = 0.0f;
    }

    for (i = 0; i < ofdm_rowsperframe; i++) {
        for (j = 0; j < ofdm_nc; j++) {
            ofdm->aphase_est_pilot_log[ofdm_nc * i + j] = 0.0f;
            ofdm->rx_amp[ofdm_nc * i + j] = 0.0f;
        }
    }

    /* default settings of options and states */

    ofdm->verbose = 0;
    ofdm->timing_en = true;
    ofdm->foff_est_en = true;
    ofdm->phase_est_en = true;
    ofdm->phase_est_bandwidth = high_bw;
    ofdm->phase_est_bandwidth_mode = AUTO_PHASE_EST;

    ofdm->foff_est_gain = 0.1f;
    ofdm->foff_est_hz = 0.0f;
    ofdm->sample_point = 0;
    ofdm->timing_est = 0;
    ofdm->timing_valid = 0;
    ofdm->timing_mx = 0.0f;
    ofdm->nin = ofdm_samplesperframe;
    ofdm->mean_amp = 0.0f;
    ofdm->foff_metric = 0.0f;

    /*
     * Unique Word symbol placement, designed to get no false syncs at any
     * freq offset.  Use ofdm_dev.m, debug_false_sync() to test.  Note we
     * need to pair the UW bits so they fit into symbols.  The LDPC decoder
     * works on symbols so we can't break up any symbols into UW/LDPC bits.
     */
    uw_ind = MALLOC(sizeof (int) * ofdm_nuwbits);
    assert(uw_ind != NULL);

    uw_ind_sym = MALLOC(sizeof (int) * (ofdm_nuwbits / 2));
    assert(uw_ind_sym != NULL);

    /*
     * The Unique Word is placed in different indexes based on
     * the number of carriers requested.
     */
    for (i = 0, j = 0; i < (ofdm_nuwbits / 2); i++, j += 2) {
        int val = floorf((i + 1) * (ofdm_nc + 1) / 2);
        uw_ind_sym[i] = val;             // symbol index

        uw_ind[j    ] = (val * 2);       // bit index 1
        uw_ind[j + 1] = (val * 2) + 1;   // bit index 2
    }

    tx_uw_syms = MALLOC(sizeof (complex float) * (ofdm_nuwbits / 2));
    assert(tx_uw_syms != NULL);

    for (i = 0; i < (ofdm_nuwbits / 2); i++) {
        tx_uw_syms[i] = 1.0f;      // qpsk_mod(0:0)
    }

    /* sync state machine */

    ofdm->sync_state = search;
    ofdm->last_sync_state = search;
    ofdm->sync_state_interleaver = search;
    ofdm->last_sync_state_interleaver = search;

    ofdm->uw_errors = 0;
    ofdm->sync_counter = 0;
    ofdm->frame_count = 0;
    ofdm->frame_count_interleaver = 0;
    ofdm->sync_start = false;
    ofdm->sync_end = false;
    ofdm->sync_mode = autosync;

    /* create the OFDM pilot time-domain waveform */

    complex float *temp = MALLOC(sizeof (complex float) * ofdm_m);
    assert(temp != NULL);

    idft(ofdm, temp, ofdm->pilots);

    /*
     * pilot_samples is 160 samples, but timing and freq offset est
     * were found by experiment to work better without a cyclic
     * prefix, so we uses zeroes instead.
     */

    /* zero out Cyclic Prefix (CP) time-domain values */

    for (i = 0; i < ofdm_ncp; i++) {
        ofdm->pilot_samples[i] = 0.0f;
    }

    /* Now copy the whole thing after the above */

    for (i = ofdm_ncp, j = 0; j < ofdm_m; i++, j++) {
        ofdm->pilot_samples[i] = temp[j];
    }

    FREE(temp);    /* finished with temp */

    /* calculate constant used to normalise timing correlation maximum */

    float acc = 0.0f;

    for (i = 0; i < (ofdm_m + ofdm_ncp); i++) {
        acc += cnormf(ofdm->pilot_samples[i]);
    }

    ofdm->timing_norm = (ofdm_m + ofdm_ncp) * acc;
    ofdm->clock_offset_counter = 0;
    ofdm->sig_var = ofdm->noise_var = 1.0f;
    ofdm->tx_bpf_en = false;
    ofdm->dpsk = false;
    
    return ofdm; /* Success */
}

static void allocate_tx_bpf(struct OFDM *ofdm) {
    //fprintf(stderr, "allocate_tx_bpf()\n");
    ofdm->ofdm_tx_bpf = MALLOC(sizeof(struct quisk_cfFilter));
    assert(ofdm->ofdm_tx_bpf != NULL);
    
    /* Transmit bandpass filter; complex coefficients, center frequency */

    quisk_filt_cfInit(ofdm->ofdm_tx_bpf, filtP550S750, sizeof (filtP550S750) / sizeof (float));
    quisk_cfTune(ofdm->ofdm_tx_bpf, ofdm_tx_centre / ofdm_fs);
}

static void deallocate_tx_bpf(struct OFDM *ofdm) {
    //fprintf(stderr, "deallocate_tx_bpf()\n");
    assert(ofdm->ofdm_tx_bpf != NULL);
    quisk_filt_destroy(ofdm->ofdm_tx_bpf);
    FREE(ofdm->ofdm_tx_bpf);
    ofdm->ofdm_tx_bpf = NULL;
}

void ofdm_destroy(struct OFDM *ofdm) {
    int i;

    if (ofdm->ofdm_tx_bpf) {
        deallocate_tx_bpf(ofdm);
    }

    FREE(ofdm->pilot_samples);
    FREE(ofdm->rxbuf);
    FREE(ofdm->pilots);

    for (i = 0; i < (ofdm_ns + 3); i++) { /* 2D array */
        FREE(ofdm->rx_sym[i]);
    }

    FREE(ofdm->rx_sym);
    FREE(ofdm->rx_np);
    FREE(ofdm->rx_amp);
    FREE(ofdm->aphase_est_pilot_log);
    FREE(ofdm->tx_uw);
    FREE(tx_uw_syms);
    FREE(uw_ind);
    FREE(uw_ind_sym);
    FREE(ofdm);
}

/*
 * Convert frequency domain into time domain
 *
 * This algorithm was designed for speed
 */
static void idft(struct OFDM *ofdm, complex float *result, complex float *vector) {
    int row, col;

    result[0] = 0.0f;

    for (col = 0; col < (ofdm_nc + 2); col++) {
        result[0] += vector[col];    // cexp(j0) == 1
    }

    result[0] *= ofdm_inv_m;

    for (row = 1; row < ofdm_m; row++) {
        complex float c = cmplx(ofdm_tx_nlower * ofdm_doc *row);
        complex float delta = cmplx(ofdm_doc * row);

        result[row] = 0.0f;

        for (col = 0; col < (ofdm_nc + 2); col++) {
            result[row] += (vector[col] * c);
            c *= delta;
        }

        result[row] *= ofdm_inv_m;
    }
}

/*
 * Convert time domain into frequency domain
 *
 * This algorithm was designed for speed
 */
static void dft(struct OFDM *ofdm, complex float *result, complex float *vector) {
    int row, col;

    for (col = 0; col < (ofdm_nc + 2); col++) {
        result[col] = vector[0];                 // conj(cexp(j0)) == 1
    }

    for (col = 0; col < (ofdm_nc + 2); col++) {
        float tval = (ofdm_rx_nlower + col) * ofdm_doc;
        complex float c = cmplxconj(tval);
        complex float delta = c;

        for (row = 1; row < ofdm_m; row++) {
            result[col] += (vector[row] * c);
            c *= delta;
        }
    }
}

static complex float vector_sum(complex float *a, int num_elements) {
    complex float sum = 0.0f;
    int i;

    for (i = 0; i < num_elements; i++) {
        sum += a[i];
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
 * Breaks when freq offset approaches +/- symbol rate (e.g
 * +/- 25 Hz for 700D).
 */
static int est_timing(struct OFDM *ofdm, complex float *rx, int length,
  int fcoarse, float *timing_mx, int *timing_valid, int step) {
    complex float corr_st, corr_en;
    int Ncorr = length - (ofdm_samplesperframe + (ofdm_m + ofdm_ncp));
    float corr[Ncorr];
    int i, j;
    float acc = 0.0f;

    for (i = 0; i < length; i++) {
        acc += cnormf(rx[i]);
    }

    float av_level = 1.0f/(2.0f * sqrtf(ofdm->timing_norm * acc / length) + 1E-12f);

    /* precompute the freq shift mulyiplied by pilot samples ouside of main loop */

    PROFILE_VAR(wvecpilot);
    PROFILE_SAMPLE(wvecpilot);

    complex float wvec_pilot[ofdm_m + ofdm_ncp];

    switch(fcoarse) {
    case -40:
      for (j = 0; j < (ofdm_m + ofdm_ncp); j++)
	wvec_pilot[j] = conjf(ofdm_wval[j]*ofdm->pilot_samples[j]);
      break;
    case 0:
      for (j = 0; j < (ofdm_m + ofdm_ncp); j++)
	wvec_pilot[j] = conjf(ofdm->pilot_samples[j]);
      break;
    case 40:
      for (j = 0; j < (ofdm_m + ofdm_ncp); j++)
	wvec_pilot[j] = ofdm_wval[j]*conjf(ofdm->pilot_samples[j]);
      break;
    default:
      assert(0);
    }

    /* use of __REAL__ provides a speed in increase of 10ms/frame during acquisition, however complex
       is fast enough for real time opration */
    
#if defined(__EMBEDDED__) && defined(__REAL__)
    float rx_real[length];
    float wvec_pilot_real[ofdm_m + ofdm_ncp];
    float wvec_pilot_imag[ofdm_m + ofdm_ncp];

    for (i = 0; i < length; i++) {
        rx_real[i] = crealf(rx[i]);
    }

    for (i = 0; i < (ofdm_m + ofdm_ncp); i++) {
        wvec_pilot_real[i] = crealf(wvec_pilot[i]);
        wvec_pilot_imag[i] = cimagf(wvec_pilot[i]);
    }
    
#endif
    PROFILE_SAMPLE_AND_LOG2(wvecpilot, "  wvecpilot");
    PROFILE_VAR(corr_start);
    PROFILE_SAMPLE(corr_start);

    for (i = 0; i < Ncorr; i += step) {
        corr_st = 0.0f;
        corr_en = 0.0f;

#ifdef __EMBEDDED__
#ifdef __REAL__
	float re,im;

	arm_dot_prod_f32(&rx_real[i], wvec_pilot_real, ofdm_m + ofdm_ncp, &re);
	arm_dot_prod_f32(&rx_real[i], wvec_pilot_imag, ofdm_m + ofdm_ncp, &im);
	corr_st = re + im * I;

	arm_dot_prod_f32(&rx_real[i+ ofdm_samplesperframe], wvec_pilot_real, ofdm_m + ofdm_ncp, &re);
	arm_dot_prod_f32(&rx_real[i+ ofdm_samplesperframe], wvec_pilot_imag, ofdm_m + ofdm_ncp, &im);
	corr_en = re + im * I;
#else
	float re,im;

	arm_cmplx_dot_prod_f32(&rx[i], wvec_pilot, ofdm_m + ofdm_ncp, &re, &im);
	corr_st = re + im * I;

	arm_cmplx_dot_prod_f32(&rx[i+ ofdm_samplesperframe], wvec_pilot, ofdm_m + ofdm_ncp, &re, &im);
	corr_en = re + im * I;
#endif        
#else	
	for (j = 0; j < (ofdm_m + ofdm_ncp); j++) {
            int ind = i + j;

	    corr_st = corr_st + (rx[ind                       ] * wvec_pilot[j]);
            corr_en = corr_en + (rx[ind + ofdm_samplesperframe] * wvec_pilot[j]);
        }
#endif	
        corr[i] = (cabsf(corr_st) + cabsf(corr_en)) * av_level;
    }

    PROFILE_SAMPLE_AND_LOG2(corr_start, "  corr");

    /* find the max magnitude and its index */

    int timing_est = 0;
    *timing_mx = 0.0f;
    
    for (i = 0; i < Ncorr; i+=step) {
        if (corr[i] > *timing_mx) {
            *timing_mx = corr[i];
            timing_est = i;
        }
    }

    // only declare timing valid if there are enough samples in rxbuf to demodulate a frame
    *timing_valid = (cabsf(rx[timing_est]) > 0.0) && (*timing_mx > ofdm_timing_mx_thresh); 
    if (ofdm->verbose > 2) {
        fprintf(stderr, "  av_level: %f  max: %f timing_est: %d timing_valid: %d\n", (double) av_level,
             (double) *timing_mx, timing_est, *timing_valid);
    }

    return timing_est;
}

/*
 * Determines frequency offset at current timing estimate, used for
 * coarse freq offset estimation during acquisition.  Works up to +/-
 * the symbol rate, e.g. +/- 25Hz for the FreeDV 700D configuration.
 */
static float est_freq_offset_pilot_corr(struct OFDM *ofdm, complex float *rx, int timing_est, int fcoarse) {
    int st = -20; int en = 20; float foff_est = 0.0f; float Cabs_max = 0.0f;

    /* precompute the freq shift mulyiplied by pilot samples ouside of main loop */

    complex float wvec_pilot[ofdm_m + ofdm_ncp];
    int j;

    switch(fcoarse) {
    case -40:
      for (j = 0; j < (ofdm_m + ofdm_ncp); j++)
	wvec_pilot[j] = conjf(ofdm_wval[j]*ofdm->pilot_samples[j]);
      break;
    case 0:
      for (j = 0; j < (ofdm_m + ofdm_ncp); j++)
	wvec_pilot[j] = conjf(ofdm->pilot_samples[j]);
      break;
    case 40:
      for (j = 0; j < (ofdm_m + ofdm_ncp); j++)
	wvec_pilot[j] = ofdm_wval[j]*conjf(ofdm->pilot_samples[j]);
      break;
    default:
      assert(0);
    }

    // sample sum of DFT magnitude of correlated signals at each freq offset and look for peak
    for (int f = st; f < en; f++) {
        complex float corr_st = 0.0f;
        complex float corr_en = 0.0f;
        float tmp = TAU * f / ofdm_fs;
	complex float delta = cmplxconj(tmp);
	complex float w = cmplxconj(0.0f);
	int i;

        for (i = 0; i < (ofdm_m + ofdm_ncp); i++) {
            // "mix" down (correlate) the pilot sequences from frame with 0 Hz offset pilot samples
            complex float csam = wvec_pilot[i] * w;
            int est = timing_est + i;

            corr_st += rx[est                       ] * csam;
            corr_en += rx[est + ofdm_samplesperframe] * csam;
	    w = w * delta;
	}

	float Cabs = cabsf(corr_st) + cabsf(corr_en);

	if (Cabs > Cabs_max) {
	  Cabs_max = Cabs;
	  foff_est = f;
	}
    }

    ofdm->foff_metric = 0.0f; // not used in this version of freq est algorithm

    if (ofdm->verbose > 2) {
        fprintf(stderr, "cabs_max: %f  foff_est: %f\n", (double) Cabs_max, (double) foff_est);
    }

    return foff_est;
}

/*
 * ----------------------------------------------
 * ofdm_txframe - modulates one frame of symbols
 * ----------------------------------------------
 */
void ofdm_txframe(struct OFDM *ofdm, complex float *tx, complex float *tx_sym_lin) {
    complex float aframe[ofdm_ns][ofdm_nc + 2];
    complex float asymbol[ofdm_m];
    complex float asymbol_cp[ofdm_m + ofdm_ncp];
    int i, j, k, m;

    /* initialize aframe to complex zero */

    for (i = 0; i < ofdm_ns; i++) {
        for (j = 0; j < (ofdm_nc + 2); j++) {
            aframe[i][j] = 0.0f;
        }
    }

    /* copy in a row of complex pilots to first row */

    for (i = 0; i < (ofdm_nc + 2); i++) {
        aframe[0][i] = ofdm->pilots[i];
    }

    /*
     * Place symbols in multi-carrier frame with pilots
     * This will place boundary values of complex zero around data
     */
    for (i = 1; i <= ofdm_rowsperframe; i++) {

        /* copy in the Nc complex values with [0 Nc 0] or (Nc + 2) total */

        for (j = 1; j < (ofdm_nc + 1); j++) {
            aframe[i][j] = tx_sym_lin[((i - 1) * ofdm_nc) + (j - 1)];
            if (ofdm->dpsk == true) {
                aframe[i][j] *= aframe[i-1][j];
            }
        }
    }

    /* OFDM up-convert symbol by symbol so we can add CP */

    for (i = 0, m = 0; i < ofdm_ns; i++, m += (ofdm_m + ofdm_ncp)) {
        idft(ofdm, asymbol, aframe[i]);

        /* Copy the last Ncp samples to the front */

        for (j = (ofdm_m - ofdm_ncp), k = 0; j < ofdm_m; j++, k++) {
            asymbol_cp[k] = asymbol[j];
        }

        /* Now copy the all samples for this row after it */

        for (j = ofdm_ncp, k = 0; k < ofdm_m; j++, k++) {
            asymbol_cp[j] = asymbol[k];
        }

        /* Now move row to the tx output */

        for (j = 0; j < (ofdm_m + ofdm_ncp); j++) {
            tx[m + j] = asymbol_cp[j];
        }
    }

    /* optional Tx Band Pass Filter */

    if (ofdm->tx_bpf_en == true) {
        assert(ofdm->ofdm_tx_bpf != NULL);
        complex float tx_filt[ofdm_samplesperframe];

        quisk_ccfFilter(tx, tx_filt, ofdm_samplesperframe, ofdm->ofdm_tx_bpf);
        memmove(tx, tx_filt, ofdm_samplesperframe * sizeof (complex float));
    }
}

struct OFDM_CONFIG *ofdm_get_config_param() {
    return &ofdm_config;
}

int ofdm_get_nin(struct OFDM *ofdm) {
    return ofdm->nin;
}

int ofdm_get_samples_per_frame() {
    return ofdm_samplesperframe;
}

int ofdm_get_max_samples_per_frame() {
    return ofdm_max_samplesperframe;
}

int ofdm_get_bits_per_frame() {
    return ofdm_bitsperframe;
}

void ofdm_set_verbose(struct OFDM *ofdm, int level) {
    ofdm->verbose = level;
}

void ofdm_set_timing_enable(struct OFDM *ofdm, bool val) {
    ofdm->timing_en = val;

    if (ofdm->timing_en == false) {
        /* manually set ideal timing instant */

        ofdm->sample_point = (ofdm_ncp - 1);
    }
}

int ofdm_get_phase_est_bandwidth_mode(struct OFDM *ofdm) {
    return ofdm->phase_est_bandwidth_mode;    /* int version of enum */
}

void ofdm_set_phase_est_bandwidth_mode(struct OFDM *ofdm, int val) {
    assert((val == AUTO_PHASE_EST) || (val == LOCKED_PHASE_EST));
    ofdm->phase_est_bandwidth_mode = val;
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
    //fprintf(stderr, "ofdm_set_tx_bpf: val: %d ofdm_tx_bpf: %p \n", val, ofdm->ofdm_tx_bpf);
    if (val == true) {
    	 if (ofdm->ofdm_tx_bpf == NULL)
             allocate_tx_bpf(ofdm);
    	ofdm->tx_bpf_en = true;
    }
    else {
    	if (ofdm->ofdm_tx_bpf != NULL)
            deallocate_tx_bpf(ofdm);
    	ofdm->tx_bpf_en = false;
    }
}

void ofdm_set_dpsk(struct OFDM *ofdm, bool val) {
    ofdm->dpsk = val;
}

/*
 * --------------------------------------
 * ofdm_mod - modulates one frame of bits
 * --------------------------------------
 */
void ofdm_mod(struct OFDM *ofdm, COMP *result, const int *tx_bits) {
    int length = ofdm_bitsperframe / ofdm_bps;
    complex float *tx = (complex float *) &result[0]; // complex has same memory layout
    complex float tx_sym_lin[length];
    int dibit[2];
    int s, i;

    if (ofdm_bps == 1) {
        /* Here we will have Nbitsperframe / 1 */

        for (s = 0; s < length; s++) {
            tx_sym_lin[s] = (float) (2 * tx_bits[s] - 1);
        }
    } else if (ofdm_bps == 2) {
        /* Here we will have Nbitsperframe / 2 */

        for (s = 0, i = 0; i < length; s += 2, i++) {
            dibit[0] = tx_bits[s + 1] & 0x1;
            dibit[1] = tx_bits[s    ] & 0x1;

            tx_sym_lin[i] = qpsk_mod(dibit);
        }
    }

    ofdm_txframe(ofdm, tx, tx_sym_lin);
}

/*
 * ----------------------------------------------------------------------------------
 * ofdm_sync_search - attempts to find coarse sync parameters for modem initial sync
 * ----------------------------------------------------------------------------------
 */

/*
 * This is a wrapper to maintain the older functionality
 * with an array of COMPs as input
 */
int ofdm_sync_search(struct OFDM *ofdm, COMP *rxbuf_in) {
    /*
     * insert latest input samples into rxbuf
     * so it is primed for when we have to call ofdm_demod()
     */

    /* note can't use memcpy when src and dest overlap */
    memmove(&ofdm->rxbuf[0], &ofdm->rxbuf[ofdm->nin],
           (ofdm_rxbuf - ofdm->nin) * sizeof (complex float));
    memmove(&ofdm->rxbuf[(ofdm_rxbuf - ofdm->nin)],
        rxbuf_in, ofdm->nin * sizeof (complex float));

    return(ofdm_sync_search_core(ofdm));
}

/*
 * This is a wrapper with a new interface to reduce memory allocated.
 * This works with ofdm_demod and freedv_api. Gain is not used here.
 */
int ofdm_sync_search_shorts(struct OFDM *ofdm, short *rxbuf_in, float gain) {
    int i, j;

    /* shift the buffer left based on nin */

    memmove(&ofdm->rxbuf[0], &ofdm->rxbuf[ofdm->nin],
            (ofdm_rxbuf - ofdm->nin) * sizeof (complex float));

    /* insert latest input samples onto tail of rxbuf */

    for (j = 0, i = (ofdm_rxbuf - ofdm->nin); i < ofdm_rxbuf; j++, i++) {
        ofdm->rxbuf[i] = ((float)rxbuf_in[j] / 32767.0f);
    }

    return ofdm_sync_search_core(ofdm);
}

/*
 * Attempts to find coarse sync parameters for modem initial sync
 */
static int ofdm_sync_search_core(struct OFDM *ofdm) {
    int act_est, afcoarse;

    /* Attempt coarse timing estimate (i.e. detect start of frame) at a range of frequency offsets */

    int st = ofdm_m + ofdm_ncp + ofdm_samplesperframe;
    int en = st + 2 * ofdm_samplesperframe + ofdm_m + ofdm_ncp;

    int fcoarse = 0;
    float atiming_mx, timing_mx = 0.0f;
    int ct_est = 0;     
    int atiming_valid, timing_valid = 0;

    PROFILE_VAR(timing_start);
    PROFILE_SAMPLE(timing_start);

    for (afcoarse = -40; afcoarse <= 40; afcoarse += 40) {
        act_est = est_timing(ofdm, &ofdm->rxbuf[st], (en - st), afcoarse, &atiming_mx, &atiming_valid, 2);

        if (atiming_mx > timing_mx) {
            ct_est = act_est;
            timing_mx = atiming_mx;
            fcoarse = afcoarse;
            timing_valid = atiming_valid;
        }
    }

    PROFILE_SAMPLE_AND_LOG2(timing_start, "  timing");

    /* refine freq est within -/+ 20 Hz window */

    PROFILE_VAR(freq_start);
    PROFILE_SAMPLE(freq_start);

    ofdm->coarse_foff_est_hz = est_freq_offset_pilot_corr(ofdm, &ofdm->rxbuf[st], ct_est, fcoarse);
    ofdm->coarse_foff_est_hz += fcoarse;

    PROFILE_SAMPLE_AND_LOG2(freq_start, "  freq");

    if (ofdm->verbose != 0) {
        fprintf(stderr, "   ct_est: %4d foff_est: %4.1f timing_valid: %d timing_mx: %5.4f\n",
                ct_est, (double) ofdm->coarse_foff_est_hz, timing_valid,
                (double)timing_mx);
    }

    if (timing_valid != 0) {
        /* potential candidate found .... */

        /* calculate number of samples we need on next buffer to get into sync */

        ofdm->nin = ct_est;

        /* reset modem states */

        ofdm->sample_point = ofdm->timing_est = 0;
        ofdm->foff_est_hz = ofdm->coarse_foff_est_hz;
        ofdm->timing_valid = timing_valid;
        ofdm->timing_mx = timing_mx;
    } else {
        ofdm->nin = ofdm_samplesperframe;
    }

    ofdm->timing_mx = timing_mx;

    return ofdm->timing_valid;
}

/*
 * ------------------------------------------
 * ofdm_demod - Demodulates one frame of bits
 * ------------------------------------------
 */

/*
 * This is a wrapper to maintain the older functionality with an
 * array of COMPs as input
 */
void ofdm_demod(struct OFDM *ofdm, int *rx_bits, COMP *rxbuf_in) {
    complex float *rx = (complex float *) &rxbuf_in[0]; // complex has same memory layout
    int i, j;

    /* shift the buffer left based on nin */
    for (i = 0, j = ofdm->nin; i < (ofdm_rxbuf - ofdm->nin); i++, j++) {
        ofdm->rxbuf[i] = ofdm->rxbuf[j];
    }

    /* insert latest input samples onto tail of rxbuf */
    for (j = 0, i = (ofdm_rxbuf - ofdm->nin); i < ofdm_rxbuf; j++, i++) {
        ofdm->rxbuf[i] = rx[j];
    }

    ofdm_demod_core(ofdm, rx_bits);
}

/*
 * This is a wrapper with a new interface to reduce memory allocated.
 * This works with ofdm_demod and freedv_api. Gain is not used here.
 */
void ofdm_demod_shorts(struct OFDM *ofdm, int *rx_bits, short *rxbuf_in, float gain) {
    int i, j;

    /* shift the buffer left based on nin */

    for (i = 0, j = ofdm->nin; i < (ofdm_rxbuf - ofdm->nin); i++, j++) {
        ofdm->rxbuf[i] = ofdm->rxbuf[j];
    }

    /* insert latest input samples onto tail of rxbuf */

    for (j = 0, i = (ofdm_rxbuf - ofdm->nin); i < ofdm_rxbuf; j++, i++) {
        ofdm->rxbuf[i] = ((float)rxbuf_in[j] / 32767.0f);
    }

    ofdm_demod_core(ofdm, rx_bits);
}

/*
 * This is the rest of the function which expects that the data is
 * already in ofdm->rxbuf
 */
static void ofdm_demod_core(struct OFDM *ofdm, int *rx_bits) {
    int prev_timing_est = ofdm->timing_est;
    int i, j, k, rr, st, en;

    /*
     * get user and calculated freq offset
     */
    float woff_est = TAU * ofdm->foff_est_hz / ofdm_fs;

    /* update timing estimate ---------------------------------------------- */

    if (ofdm->timing_en == true) {
        /* update timing at start of every frame */

        st = ((ofdm_m + ofdm_ncp) + ofdm_samplesperframe) - floorf(ofdm_ftwindowwidth / 2) + ofdm->timing_est;
        en = st + ofdm_samplesperframe - 1 + (ofdm_m + ofdm_ncp) + ofdm_ftwindowwidth;

        complex float work[(en - st)];

        /*
         * Adjust for the frequency error by shifting the phase
         * using a conjugate multiply
         */
        for (j = 0, i = st; i < en; j++, i++) {
            work[j] = ofdm->rxbuf[i] * cmplxconj(woff_est * i);
        }

        int ft_est = est_timing(ofdm, work, (en - st), 0.0f, &ofdm->timing_mx, &ofdm->timing_valid, 1);
        
        ofdm->timing_est += ft_est - ceilf((float)ofdm_ftwindowwidth / 2) + 1;

        if (ofdm->verbose > 2) {
            fprintf(stderr, "  ft_est: %2d timing_est: %2d sample_point: %2d\n", ft_est, ofdm->timing_est,
                ofdm->sample_point);
        }

        /* Black magic to keep sample_point inside cyclic prefix.  Or something like that. */

        ofdm->sample_point = max(ofdm->timing_est + (ofdm_ncp / 4), ofdm->sample_point);
        ofdm->sample_point = min(ofdm->timing_est + ofdm_ncp, ofdm->sample_point);
    }

    /*
     * Convert the time-domain samples to the frequency-domain using the rx_sym
     * data matrix. This will be  Nc+2 carriers of 11 symbols.
     *
     * You will notice there are Nc+2 BPSK symbols for each pilot symbol, and
     * that there are Nc QPSK symbols for each data symbol.
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
    for (i = 0; i < (ofdm_ns + 3); i++) {
        for (j = 0; j < (ofdm_nc + 2); j++) {
            ofdm->rx_sym[i][j] = 0.0f;
        }
    }

    /*
     * "Previous" pilot symbol is one modem frame above.
     */
    st = (ofdm_m + ofdm_ncp) + 1 + ofdm->sample_point;
    en = st + ofdm_m;

    complex float work[ofdm_m];

    /* down-convert at current timing instant------------------------------- */

    for (k = 0, j = st; j < en; k++, j++) {
        work[k] = ofdm->rxbuf[j] * cmplxconj(woff_est * j);
    }

    /*
     * Each symbol is of course (ofdm_m + ofdm_ncp) samples long and
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
     * "This" pilot comes after the extra symbol allotted at the top, and after
     * the "previous" pilot and previous data symbols (let's call it, the previous
     * modem frame).
     *
     * So we will now be starting at "this" pilot symbol, and continuing to the
     * "next" pilot symbol.
     *
     * In this routine we also process the current data symbols.
     */
    for (rr = 0; rr < (ofdm_ns + 1); rr++) {
        st = (ofdm_m + ofdm_ncp) + ofdm_samplesperframe + (rr * (ofdm_m + ofdm_ncp)) + 1 + ofdm->sample_point;
        en = st + ofdm_m;

        /* down-convert at current timing instant---------------------------------- */

        for (k = 0, j = st; j < en; k++, j++) {
            work[k] = ofdm->rxbuf[j] * cmplxconj(woff_est * j);
        }

        /*
         * We put these Nc+2 carrier symbols into our matrix after the previous pilot:
         *
         * 1 .................. Nc+2
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
    st = (ofdm_m + ofdm_ncp) + (3 * ofdm_samplesperframe) + 1 + ofdm->sample_point;
    en = st + ofdm_m;

    /* down-convert at current timing instant------------------------------- */

    for (k = 0, j = st; j < en; k++, j++) {
        work[k] = ofdm->rxbuf[j] * cmplxconj(woff_est * j);
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
    dft(ofdm, ofdm->rx_sym[ofdm_ns + 2], work);

    /*
     * We are finished now with the DFT and down conversion
     * From here on down we are in the frequency domain
     */

    /* est freq err based on all carriers ---------------------------------- */

    if (ofdm->foff_est_en == true) {
        /*
         * sym[1] is 'this' pilot symbol, sym[9] is 'next' pilot symbol.
         *
         * By subtracting the two averages of these pilots, we find the frequency
         * by the change in phase over time.
         */
        complex float freq_err_rect =
                conjf(vector_sum(ofdm->rx_sym[1], ofdm_nc + 2)) *
                vector_sum(ofdm->rx_sym[ofdm_ns + 1], ofdm_nc + 2);

        /* prevent instability in atan(im/re) when real part near 0 */

        freq_err_rect += 1E-6f;

        float freq_err_hz = cargf(freq_err_rect) * ofdm_rs / (TAU * ofdm_ns);
        ofdm->foff_est_hz += (ofdm->foff_est_gain * freq_err_hz);
    }

    /* OK - now estimate and correct pilot phase  -------------------------- */

    complex float aphase_est_pilot_rect;
    float aphase_est_pilot[ofdm_nc + 2];
    float aamp_est_pilot[ofdm_nc + 2];

    for (i = 0; i < (ofdm_nc + 2); i++) {
        aphase_est_pilot[i] = 10.0f;
        aamp_est_pilot[i] = 0.0f;
    }

    for (i = 1; i < (ofdm_nc + 1); i++) { /* ignore first and last carrier for count */
        if (ofdm->phase_est_bandwidth == low_bw) {
            complex float symbol[3];

            /*
             * Use all pilots normally, results in low SNR performance,
             * but will fall over in high Doppler propagation
             *
             * Basically we divide the Nc+2 pilots into groups of 3
             * Then average the phase surrounding each of the data symbols.
             */
            for (k = 0, j = (i - 1); k < 3; k++, j++) {
                symbol[k] = ofdm->rx_sym[1][j] * conjf(ofdm->pilots[j]); /* this pilot conjugate */
            }

            aphase_est_pilot_rect = vector_sum(symbol, 3);

            for (k = 0, j = (i - 1); k < 3; k++, j++) {
                symbol[k] = ofdm->rx_sym[ofdm_ns + 1][j] * conjf(ofdm->pilots[j]); /* next pilot conjugate */
            }

            aphase_est_pilot_rect += vector_sum(symbol, 3);

            /* use pilots in past and future */

            for (k = 0, j = (i - 1); k < 3; k++, j++) {
                symbol[k] = ofdm->rx_sym[0][j] * conjf(ofdm->pilots[j]); /* previous pilot */
            }

            aphase_est_pilot_rect += vector_sum(symbol, 3);

            for (k = 0, j = (i - 1); k < 3; k++, j++) {
                symbol[k] = ofdm->rx_sym[ofdm_ns + 2][j] * conjf(ofdm->pilots[j]); /* future pilot */
            }

            aphase_est_pilot_rect += vector_sum(symbol, 3);
            aphase_est_pilot[i] = cargf(aphase_est_pilot_rect);

            /* amplitude is estimated over 12 pilots */

            aamp_est_pilot[i] = cabsf(aphase_est_pilot_rect / 12.0f);
        } else {
            assert(ofdm->phase_est_bandwidth == high_bw);

            /*
             * Use only symbols at 'this' and 'next' to quickly track changes
             * in phase due to high Doppler spread in propagation (no neighbor averaging).
             *
             * As less pilots are averaged, low SNR performance will be poorer
             */
            aphase_est_pilot_rect = ofdm->rx_sym[1][i] * conjf(ofdm->pilots[i]);            /* "this" pilot conjugate */
            aphase_est_pilot_rect += ofdm->rx_sym[ofdm_ns + 1][i] * conjf(ofdm->pilots[i]); /* "next" pilot conjugate */
            aphase_est_pilot[i] = cargf(aphase_est_pilot_rect);

            /* amplitude is estimated over 2 pilots */

            aamp_est_pilot[i] = cabsf(aphase_est_pilot_rect / 2.0f);
        }
    }

    /*
     * correct the phase offset using phase estimate, and demodulate
     * bits, separate loop as it runs across cols (carriers) to get
     * frame bit ordering correct
     */
    complex float rx_corr;
    int abit[2];
    int bit_index = 0;
    float sum_amp = 0.0f;

    for (rr = 0; rr < ofdm_rowsperframe; rr++) {
        /*
         * Note the i starts with the second carrier, ends with Nc+1.
         * so we ignore the first and last carriers.
         *
         * Also note we are using sym[2..8] or the seven data symbols.
         */
        for (i = 1; i < (ofdm_nc + 1); i++) {
            if (ofdm->phase_est_en == true) {
                if (ofdm->dpsk == true) {
                    /* differential detection, using pilot as reference at start of frame */
                    rx_corr = ofdm->rx_sym[rr + 2][i] * cmplxconj(cargf(ofdm->rx_sym[rr + 1][i]));
                } else  {
                    /* regular coherent detection */
                    rx_corr = ofdm->rx_sym[rr + 2][i] * cmplxconj(aphase_est_pilot[i]);
                }
            } else {
                rx_corr = ofdm->rx_sym[rr + 2][i];
            }

            /*
             * Output complex data symbols after phase correction;
             * (_np = no pilots) the pilot symbols have been removed
             */
            ofdm->rx_np[(rr * ofdm_nc) + (i - 1)] = rx_corr;

            /*
             * Note even though amp ests are the same for each col,
             * the FEC decoder likes to have one amplitude per symbol
             * so convenient to log them all
             */
            ofdm->rx_amp[(rr * ofdm_nc) + (i - 1)] = aamp_est_pilot[i];
            sum_amp += aamp_est_pilot[i];

            /*
             * Note like amps in this implementation phase ests are the
             * same for each col, but we log them for each symbol anyway
             */
            ofdm->aphase_est_pilot_log[(rr * ofdm_nc) + (i - 1)] = aphase_est_pilot[i];

            if (ofdm_bps == 1) {
                rx_bits[bit_index++] = crealf(rx_corr) > 0.0f;
            } else if (ofdm_bps == 2) {
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

    ofdm->mean_amp = 0.9f * ofdm->mean_amp + 0.1f * sum_amp / (ofdm_rowsperframe * ofdm_nc);

    /* Adjust nin to take care of sample clock offset */

    ofdm->nin = ofdm_samplesperframe;

    if (ofdm->timing_en == true) {
        ofdm->clock_offset_counter += (prev_timing_est - ofdm->timing_est);

        int thresh = (ofdm_m + ofdm_ncp) / 8;
        int tshift = (ofdm_m + ofdm_ncp) / 4;

        if (ofdm->timing_est > thresh) {
            ofdm->nin = ofdm_samplesperframe + tshift;
            ofdm->timing_est -= tshift;
            ofdm->sample_point -= tshift;
        } else if (ofdm->timing_est < -thresh) {
            ofdm->nin = ofdm_samplesperframe - tshift;
            ofdm->timing_est += tshift;
            ofdm->sample_point += tshift;
        }
    }

    /*
     * estimate signal and noise power, see ofdm_lib.m,
     * cohpsk.m for more info
     */
    complex float *rx_np = ofdm->rx_np;

    float sig_var = 0.0f;

    /*
     * sig_var gets a little large, so tamp it down
     * each step
     */
    float step = (1.0f / (ofdm_rowsperframe * ofdm_nc));

    for (i = 0; i < (ofdm_rowsperframe * ofdm_nc); i++) {
        sig_var += (cnormf(rx_np[i]) * step);
    }

    float sig_rms = sqrtf(sig_var);

    float sum_x = 0.0f;
    float sum_xx = 0.0f;
    int n = 0;

    for (i = 0; i < (ofdm_rowsperframe * ofdm_nc); i++) {
        complex float s = rx_np[i];

        if (fabsf(crealf(s)) > sig_rms) {
            sum_x += cimagf(s);
            sum_xx += cimagf(s) * cimagf(s);
            n++;
        }
    }

    /*
     * with large interfering carriers this alg can break down - in
     * that case set a benign value for noise_var that will produce a
     * sensible (probably low) SNR est
     */
    float noise_var = 1.0f;

    if (n > 1) {
        noise_var = (n * sum_xx - sum_x * sum_x) / (n * (n - 1));
    }

    ofdm->noise_var = 2.0f * noise_var;
    ofdm->sig_var = sig_var;
}

/*
 * iterate state machine
 */
void ofdm_sync_state_machine(struct OFDM *ofdm, uint8_t *rx_uw) {
    int i;

    State next_state = ofdm->sync_state;

    ofdm->sync_start = false;
    ofdm->sync_end = false;

    if (ofdm->sync_state == search) {
        if (ofdm->timing_valid) {
            ofdm->frame_count = 0;
            ofdm->sync_counter = 0;
            ofdm->sync_start = true;
            ofdm->clock_offset_counter = 0;
            next_state = trial;
        }
    }

    if ((ofdm->sync_state == synced) || (ofdm->sync_state == trial)) {
        ofdm->frame_count++;
        ofdm->frame_count_interleaver++;

        /*
         * freq offset est may be too far out, and has aliases every 1/Ts, so
         * we use a Unique Word to get a really solid indication of sync.
         */
        ofdm->uw_errors = 0;

        for (i = 0; i < ofdm_nuwbits; i++) {
            ofdm->uw_errors += ofdm->tx_uw[i] ^ rx_uw[i];
        }

        /*
         * during trial sync we don't tolerate errors so much, we look
         * for 3 consecutive frames with low error rate to confirm sync
         */
        if (ofdm->sync_state == trial) {
            if (ofdm->uw_errors > 2) {
                /* if we exceed thresh stay in trial sync */

                ofdm->sync_counter++;
                ofdm->frame_count = 0;
            }

            if (ofdm->sync_counter == 2) {
                /* if we get two bad frames drop sync and start again */

                next_state = search;
                ofdm->sync_state_interleaver = search;
                ofdm->phase_est_bandwidth = high_bw;
            }

            if (ofdm->frame_count == 4) {
                /* three good frames, sync is OK! */

                next_state = synced;
                /* change to low bandwidth, but more accurate phase estimation */
                /* but only if not locked to high */

                if (ofdm->phase_est_bandwidth_mode != LOCKED_PHASE_EST) {
                    ofdm->phase_est_bandwidth = low_bw;
                }
            }
        }

        /* once we have synced up we tolerate a higher error rate to wait out fades */

        if (ofdm->sync_state == synced) {
            if (ofdm->uw_errors > 2) {
                ofdm->sync_counter++;
            } else {
                ofdm->sync_counter = 0;
            }

            if ((ofdm->sync_mode == autosync) && (ofdm->sync_counter > 6)) {
                /* run of consecutive bad frames ... drop sync */

                next_state = search;
                ofdm->sync_state_interleaver = search;
                ofdm->phase_est_bandwidth = high_bw;
            }
        }
    }

    ofdm->last_sync_state = ofdm->sync_state;
    ofdm->last_sync_state_interleaver = ofdm->sync_state_interleaver;
    ofdm->sync_state = next_state;
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
    assert(ofdm != NULL);

    switch (sync_cmd) {
        case UN_SYNC:
            /*
             * force manual unsync, in case operator detects false sync,
             * which will cause sync state machine to have another go at sync
             */
            ofdm->sync_state = search;
            ofdm->sync_state_interleaver = search;
            break;
        case AUTO_SYNC:
            /* normal operating mode - sync state machine decides when to unsync */

            ofdm->sync_mode = autosync;
            break;
        case MANUAL_SYNC:
            /*
             * allow sync state machine to sync, but not to unsync, the
             * operator will decide that manually
             */
            ofdm->sync_mode = manualsync;
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

void ofdm_get_demod_stats(struct OFDM *ofdm, struct MODEM_STATS *stats) {
    stats->Nc = ofdm_nc;
    assert(stats->Nc <= MODEM_STATS_NC_MAX);

    float snr_est = 10.0f * log10f((0.1f + (ofdm->sig_var / ofdm->noise_var)) *
            ofdm_nc * ofdm_rs / 3000.0f);
    float total = ofdm->frame_count * ofdm_samplesperframe;

    /* fast attack, slow decay */
    if (snr_est > stats->snr_est)
        stats->snr_est = snr_est;
    else
        stats->snr_est = 0.9f * stats->snr_est + 0.1f * snr_est;
    stats->sync = ((ofdm->sync_state == synced) || (ofdm->sync_state == trial));
    stats->foff = ofdm->foff_est_hz;
    stats->rx_timing = ofdm->timing_est;
    stats->clock_offset = 0;

    if (total != 0.0f) {
        stats->clock_offset = ofdm->clock_offset_counter / total;
    }

    stats->sync_metric = ofdm->timing_mx;

#ifndef __EMBEDDED__
    assert(ofdm_rowsperframe < MODEM_STATS_NR_MAX);
    stats->nr = ofdm_rowsperframe;

    for (int c = 0; c < ofdm_nc; c++) {
        for (int r = 0; r < ofdm_rowsperframe; r++) {
            complex float rot = ofdm->rx_np[r * c] * cmplx(ROT45);

            stats->rx_symbols[r][c].real = crealf(rot);
            stats->rx_symbols[r][c].imag = cimagf(rot);
        }
    }
#endif    
}

/*
 * Assemble modem frame of bits from UW, payload bits, and txt bits
 */
void ofdm_assemble_modem_frame(struct OFDM *ofdm, uint8_t modem_frame[],
        uint8_t payload_bits[], uint8_t txt_bits[]) {
    int s, t;

    int p = 0;
    int u = 0;

    for (s = 0; s < (ofdm_bitsperframe - ofdm_ntxtbits); s++) {
        if ((u < ofdm_nuwbits) && (s == uw_ind[u])) {
            modem_frame[s] = ofdm->tx_uw[u++];
        } else {
            modem_frame[s] = payload_bits[p++];
        }
    }

    assert(u == ofdm_nuwbits);
    assert(p == (ofdm_bitsperframe - ofdm_nuwbits - ofdm_ntxtbits));

    for (t = 0; s < ofdm_bitsperframe; s++, t++) {
        modem_frame[s] = txt_bits[t];
    }

    assert(t == ofdm_ntxtbits);
}

/*
 * Assemble modem frame from UW, payload symbols, and txt bits
 */
void ofdm_assemble_modem_frame_symbols(complex float modem_frame[],
  COMP payload_syms[], uint8_t txt_bits[]) {
    complex float *payload = (complex float *) &payload_syms[0]; // complex has same memory layout
    int Nsymsperframe = ofdm_bitsperframe / ofdm_bps;
    int Nuwsyms = ofdm_nuwbits / ofdm_bps;
    int Ntxtsyms = ofdm_ntxtbits / ofdm_bps;
    int dibit[2];
    int s, t;

    int p = 0;
    int u = 0;

    for (s = 0; s < (Nsymsperframe - Ntxtsyms); s++) {
        if ((u < Nuwsyms) && (s == uw_ind_sym[u])) {
            modem_frame[s] = tx_uw_syms[u++];
        } else {
            modem_frame[s] = payload[p++];
        }
    }

    assert(u == Nuwsyms);
    assert(p == (Nsymsperframe - Nuwsyms - Ntxtsyms));

    for (t = 0; s < Nsymsperframe; s++, t += ofdm_bps) {
        dibit[1] = txt_bits[t    ] & 0x1;
        dibit[0] = txt_bits[t + 1] & 0x1;
        modem_frame[s] = qpsk_mod(dibit);
    }

    assert(t == ofdm_ntxtbits);
}

void ofdm_disassemble_modem_frame(struct OFDM *ofdm, uint8_t rx_uw[],
  COMP codeword_syms[], float codeword_amps[], short txt_bits[]) {
    complex float *codeword = (complex float *) &codeword_syms[0]; // complex has same memory layout
    int Nsymsperframe = ofdm_bitsperframe / ofdm_bps;
    int Nuwsyms = ofdm_nuwbits / ofdm_bps;
    int Ntxtsyms = ofdm_ntxtbits / ofdm_bps;
    int dibit[2];
    int s, t;

    int p = 0;
    int u = 0;

    for (s = 0; s < (Nsymsperframe - Ntxtsyms); s++) {
        if ((u < Nuwsyms) && (s == uw_ind_sym[u])) {
            qpsk_demod(ofdm->rx_np[s], dibit);

            rx_uw[ofdm_bps * u    ] = dibit[1];
            rx_uw[ofdm_bps * u + 1] = dibit[0];
            u++;
        } else {
            codeword[p] = ofdm->rx_np[s];
            codeword_amps[p] = ofdm->rx_amp[s];
            p++;
        }
    }

    assert(u == Nuwsyms);
    assert(p == (Nsymsperframe - Nuwsyms - Ntxtsyms));

    for (t = 0; s < Nsymsperframe; s++, t += ofdm_bps) {
        qpsk_demod(ofdm->rx_np[s], dibit);

        txt_bits[t    ] = dibit[1];
        txt_bits[t + 1] = dibit[0];
    }

    assert(t == ofdm_ntxtbits);
}

/*
 * Pseudo-random number generator that we can implement in C with
 * identical results to Octave.  Returns an unsigned int between 0
 * and 32767.  Used for generating test frames of various lengths.
 */
void ofdm_rand(uint16_t r[], int n) {
    uint64_t seed = 1;
    int i;

    for (i = 0; i < n; i++) {
        seed = (1103515245l * seed + 12345) % 32768;
        r[i] = seed;
    }
}

void ofdm_generate_payload_data_bits(uint8_t payload_data_bits[], int data_bits_per_frame) {
    uint16_t r[data_bits_per_frame];
    int i;

    /* construct payload data bits */

    ofdm_rand(r, data_bits_per_frame);

    for (i = 0; i < data_bits_per_frame; i++) {
        payload_data_bits[i] = r[i] > 16384;
    }
}

void ofdm_print_info(struct OFDM *ofdm) {
    char *syncmode[] = {
        "unsync",
        "autosync",
        "manualsync"
    };
    char *phase_est_bandwidth_mode[] = {
        "auto",
        "locked_high"
    };

    fprintf(stderr, "ofdm_tx_centre = %g\n", (double)ofdm_tx_centre);
    fprintf(stderr, "ofdm_rx_centre = %g\n", (double)ofdm_rx_centre);
    fprintf(stderr, "ofdm_fs = %g\n", (double)ofdm_fs);
    fprintf(stderr, "ofdm_ts = %g\n", (double)ofdm_ts);
    fprintf(stderr, "ofdm_rs = %g\n", (double)ofdm_rs);
    fprintf(stderr, "ofdm_tcp = %g\n", (double)ofdm_tcp);
    fprintf(stderr, "ofdm_inv_m = %g\n", (double)ofdm_inv_m);
    fprintf(stderr, "ofdm_tx_nlower = %g\n", (double)ofdm_tx_nlower);
    fprintf(stderr, "ofdm_rx_nlower = %g\n", (double)ofdm_rx_nlower);
    fprintf(stderr, "ofdm_doc = %g\n", (double)ofdm_doc);
    fprintf(stderr, "ofdm_timing_mx_thresh = %g\n", (double)ofdm_timing_mx_thresh);
    fprintf(stderr, "ofdm_nc = %d\n", ofdm_nc);
    fprintf(stderr, "ofdm_ns = %d\n", ofdm_ns);
    fprintf(stderr, "ofdm_bps = %d\n", ofdm_bps);
    fprintf(stderr, "ofdm_m = %d\n", ofdm_m);
    fprintf(stderr, "ofdm_ncp = %d\n", ofdm_ncp);
    fprintf(stderr, "ofdm_ftwindowwidth = %d\n", ofdm_ftwindowwidth);
    fprintf(stderr, "ofdm_bitsperframe = %d\n", ofdm_bitsperframe);
    fprintf(stderr, "ofdm_rowsperframe = %d\n", ofdm_rowsperframe);
    fprintf(stderr, "ofdm_samplesperframe = %d\n", ofdm_samplesperframe);
    fprintf(stderr, "ofdm_max_samplesperframe = %d\n", ofdm_max_samplesperframe);
    fprintf(stderr, "ofdm_rxbuf = %d\n", ofdm_rxbuf);
    fprintf(stderr, "ofdm_ntxtbits = %d\n", ofdm_ntxtbits);
    fprintf(stderr, "ofdm_nuwbits = %d\n", ofdm_nuwbits);
    fprintf(stderr, "ofdm->foff_est_gain = %g\n", (double)ofdm->foff_est_gain);
    fprintf(stderr, "ofdm->foff_est_hz = %g\n", (double)ofdm->foff_est_hz);
    fprintf(stderr, "ofdm->timing_mx = %g\n", (double)ofdm->timing_mx);
    fprintf(stderr, "ofdm->coarse_foff_est_hz = %g\n", (double)ofdm->coarse_foff_est_hz);
    fprintf(stderr, "ofdm->timing_norm = %g\n", (double)ofdm->timing_norm);
    fprintf(stderr, "ofdm->sig_var = %g\n", (double)ofdm->sig_var);
    fprintf(stderr, "ofdm->noise_var = %g\n", (double)ofdm->noise_var);
    fprintf(stderr, "ofdm->mean_amp = %g\n", (double)ofdm->mean_amp);
    fprintf(stderr, "ofdm->clock_offset_counter = %d\n", ofdm->clock_offset_counter);
    fprintf(stderr, "ofdm->verbose = %d\n", ofdm->verbose);
    fprintf(stderr, "ofdm->sample_point = %d\n", ofdm->sample_point);
    fprintf(stderr, "ofdm->timing_est = %d\n", ofdm->timing_est);
    fprintf(stderr, "ofdm->timing_valid = %d\n", ofdm->timing_valid);
    fprintf(stderr, "ofdm->nin = %d\n", ofdm->nin);
    fprintf(stderr, "ofdm->uw_errors = %d\n", ofdm->uw_errors);
    fprintf(stderr, "ofdm->sync_counter = %d\n", ofdm->sync_counter);
    fprintf(stderr, "ofdm->frame_count = %d\n", ofdm->frame_count);
    fprintf(stderr, "ofdm->sync_start = %s\n", ofdm->sync_start ? "true" : "false");
    fprintf(stderr, "ofdm->sync_end = %s\n", ofdm->sync_end ? "true" : "false");
    fprintf(stderr, "ofdm->sync_mode = %s\n", syncmode[ofdm->sync_mode]);
    fprintf(stderr, "ofdm->frame_count_interleaver = %d\n", ofdm->frame_count_interleaver);
    fprintf(stderr, "ofdm->timing_en = %s\n", ofdm->timing_en ? "true" : "false");
    fprintf(stderr, "ofdm->foff_est_en = %s\n", ofdm->foff_est_en ? "true" : "false");
    fprintf(stderr, "ofdm->phase_est_en = %s\n", ofdm->phase_est_en ? "true" : "false");
    fprintf(stderr, "ofdm->tx_bpf_en = %s\n", ofdm->tx_bpf_en ? "true" : "false");
    fprintf(stderr, "ofdm->dpsk = %s\n", ofdm->dpsk ? "true" : "false");
    fprintf(stderr, "ofdm->phase_est_bandwidth_mode = %s\n", phase_est_bandwidth_mode[ofdm->phase_est_bandwidth_mode]);
}
