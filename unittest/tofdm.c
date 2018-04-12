/*---------------------------------------------------------------------------*\

  FILE........: tofdm.c
  AUTHORS.....: David Rowe & Steve Sampson
  DATE CREATED: June 2017

  Tests for the C version of the OFDM modem.  This program
  outputs a file of Octave vectors that are loaded and automatically
  tested against the Octave version of the modem by the Octave script
  tofdm.m

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

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>

#include "ofdm_internal.h"
#include "codec2_ofdm.h"
#include "octave.h"
#include "test_bits_ofdm.h"
#include "comp_prim.h"

#define NFRAMES 30
#define SAMPLE_CLOCK_OFFSET_PPM 100
#define FOFF_HZ 0.5f

#define ASCALE  (2E5*1.1491/2.0) /* scale from shorts back to floats */

/*---------------------------------------------------------------------------*\

  FUNCTION....: fs_offset()
  AUTHOR......: David Rowe
  DATE CREATED: May 2015

  Simulates small Fs offset between mod and demod.
  (Note: Won't work with float, works OK with double)

\*---------------------------------------------------------------------------*/

static int fs_offset(COMP out[], COMP in[], int n, float sample_rate_ppm) {
    double f;
    int t1, t2;

    double tin = 0.0;
    int tout = 0;

    while (tin < (double) (n-1)) {
      t1 = (int) floor(tin);
      t2 = (int) ceil(tin);
      assert(t2 < n);

      f = (tin - (double) t1);

      out[tout].real = (1.0 - f) * in[t1].real + f * in[t2].real;
      out[tout].imag = (1.0 - f) * in[t1].imag + f * in[t2].imag;

      tout += 1;
      tin  += 1.0 + sample_rate_ppm / 1E6;
    }
    //printf("n: %d tout: %d tin: %f\n", n, tout, tin);
    
    return tout;
}

#ifndef ARM_MATH_CM4
  #define SINF(a) sinf(a)
  #define COSF(a) cosf(a)
#else
  #define SINF(a) arm_sin_f32(a)
  #define COSF(a) arm_cos_f32(a)
#endif

/*---------------------------------------------------------------------------*\

  FUNCTION....: freq_shift()
  AUTHOR......: David Rowe
  DATE CREATED: 26/4/2012

  Frequency shift modem signal.  The use of complex input and output allows
  single sided frequency shifting (no images).

\*---------------------------------------------------------------------------*/

static void freq_shift(COMP rx_fdm_fcorr[], COMP rx_fdm[], float foff, COMP *foff_phase_rect, int nin) {
    float temp = (TAU * foff / OFDM_FS);
    COMP  foff_rect = { COSF(temp), SINF(temp) };
    int   i;

    for (i = 0; i < nin; i++) {
	*foff_phase_rect = cmult(*foff_phase_rect, foff_rect);
	rx_fdm_fcorr[i] = cmult(rx_fdm[i], *foff_phase_rect);
    }

    /* normalise digital oscillator as the magnitude can drift over time */

    float mag = cabsolute(*foff_phase_rect);
    foff_phase_rect->real /= mag;
    foff_phase_rect->imag /= mag;
}

int main(int argc, char *argv[])
{
    int            samples_per_frame = ofdm_get_samples_per_frame();
    int            max_samples_per_frame = ofdm_get_max_samples_per_frame();

    struct OFDM   *ofdm;
    COMP           tx[samples_per_frame];         /* one frame of tx samples */

    int            rx_bits[OFDM_BITSPERFRAME];    /* one frame of rx bits    */

    /* log arrays */

    int            tx_bits_log[OFDM_BITSPERFRAME*NFRAMES];
    COMP           tx_log[samples_per_frame*NFRAMES];
    COMP           rx_log[samples_per_frame*NFRAMES];
    COMP           rxbuf_in_log[max_samples_per_frame*NFRAMES];
    COMP           rxbuf_log[OFDM_RXBUF*NFRAMES];
    COMP           rx_sym_log[(OFDM_NS + 3)*NFRAMES][OFDM_NC + 2];
    float          phase_est_pilot_log[OFDM_ROWSPERFRAME*NFRAMES][OFDM_NC];
    COMP           rx_np_log[OFDM_ROWSPERFRAME*OFDM_NC*NFRAMES];
    float          rx_amp_log[OFDM_ROWSPERFRAME*OFDM_NC*NFRAMES];
    float          foff_hz_log[NFRAMES];
    int            rx_bits_log[OFDM_BITSPERFRAME*NFRAMES];
    int            timing_est_log[NFRAMES];
    int            timing_valid_log[NFRAMES];
    float          timing_mx_log[NFRAMES];
    float          coarse_foff_est_hz_log[NFRAMES];
    int            sample_point_log[NFRAMES];

    FILE          *fout;
    int            f,i,j;

    ofdm = ofdm_create(OFDM_CONFIG_700D);
    assert(ofdm != NULL);
    
    /* Main Loop ---------------------------------------------------------------------*/

    for(f=0; f<NFRAMES; f++) {

	/* --------------------------------------------------------*\
	                          Mod
	\*---------------------------------------------------------*/

        /* todo: add a longer sequence of test bits through
           ofdm_get/put test bits functin similat to cohpsk/fdmdv */

        ofdm_mod(ofdm, (COMP*)tx, test_bits_ofdm);

        /* tx vector logging */

	memcpy(&tx_bits_log[OFDM_BITSPERFRAME*f], test_bits_ofdm, sizeof(int)*OFDM_BITSPERFRAME);
	memcpy(&tx_log[samples_per_frame*f], tx, sizeof(COMP)*samples_per_frame);
    }

    /* --------------------------------------------------------*\
	                        Channel
    \*---------------------------------------------------------*/

    fs_offset(rx_log, tx_log, samples_per_frame*NFRAMES, SAMPLE_CLOCK_OFFSET_PPM);

    COMP foff_phase_rect = {1.0f, 0.0f};

    freq_shift(rx_log, rx_log, FOFF_HZ, &foff_phase_rect, samples_per_frame * NFRAMES);

    /* --------------------------------------------------------*\
	                        Demod
    \*---------------------------------------------------------*/

    /* Init/pre-load rx with ideal timing so we can test with timing estimation disabled */

    int  Nsam = samples_per_frame*NFRAMES;
    int  prx = 0;
    int  nin = samples_per_frame + 2*(OFDM_M+OFDM_NCP);

    int  lnew;
    COMP rxbuf_in[max_samples_per_frame];

    #define FRONT_LOAD
    #ifdef FRONT_LOAD
    for (i=0; i<nin; i++,prx++) {
         ofdm->rxbuf[OFDM_RXBUF-nin+i] = rx_log[prx].real + I*rx_log[prx].imag;
    }
    #endif
    
    int nin_tot = 0;

    /* disable estimators for initial testing */

    ofdm_set_verbose(ofdm, false);
    ofdm_set_timing_enable(ofdm, true);
    ofdm_set_foff_est_enable(ofdm, true);
    ofdm_set_phase_est_enable(ofdm, true);

    //#define TESTING_FILE
    #ifdef TESTING_FILE
    FILE *fin=fopen("/home/david/codec2-dev/octave/ofdm_test.raw", "rb");
    assert(fin != NULL);
    int Nbitsperframe = ofdm_get_bits_per_frame(ofdm);
    int Nmaxsamperframe = ofdm_get_max_samples_per_frame();
    short rx_scaled[Nmaxsamperframe];
    #endif
    
    for(f=0; f<NFRAMES; f++) {
        /* For initial testing, timing est is off, so nin is always
           fixed.  TODO: we need a constant for rxbuf_in[] size that
           is the maximum possible nin */

        nin = ofdm_get_nin(ofdm);
        assert(nin <= max_samples_per_frame);

        /* Insert samples at end of buffer, set to zero if no samples
           available to disable phase estimation on future pilots on
           last frame of simulation. */

        if ((Nsam-prx) < nin) {
            lnew = Nsam-prx;
        } else {
            lnew = nin;
        }
        //printf("nin: %d prx: %d lnew: %d\n", nin, prx, lnew);
        for(i=0; i<nin; i++) {
            rxbuf_in[i].real = 0.0;
            rxbuf_in[i].imag = 0.0;
        }

        if (lnew) {
            for(i=0; i<lnew; i++, prx++) {
                rxbuf_in[i] = rx_log[prx];
            }
        }
        assert(prx <= max_samples_per_frame*NFRAMES);

        #ifdef TESTING_FILE
        fread(rx_scaled, sizeof(short), nin, fin);
        for(i=0; i<nin; i++) {
	    rxbuf_in[i].real = (float)rx_scaled[i]/ASCALE;
            rxbuf_in[i].imag = 0.0;
        }
        #endif

        ofdm_demod(ofdm, rx_bits, rxbuf_in);
        
        #ifdef TESTING_FILE
        int Nerrs = 0;
        for(i=0; i<Nbitsperframe; i++) {
            if (test_bits_ofdm[i] != rx_bits[i]) {
                Nerrs++;
            }
        }
        printf("f: %d Nerr: %d\n", f, Nerrs);
        #endif
        
        /* rx vector logging -----------------------------------*/

        assert(nin_tot < samples_per_frame*NFRAMES);
	memcpy(&rxbuf_in_log[nin_tot], rxbuf_in, sizeof(COMP)*nin);
        nin_tot += nin;

        for(i=0; i<OFDM_RXBUF; i++) {
            rxbuf_log[OFDM_RXBUF*f+i].real = crealf(ofdm->rxbuf[i]);
            rxbuf_log[OFDM_RXBUF*f+i].imag = cimagf(ofdm->rxbuf[i]);
        }

        for (i = 0; i < (OFDM_NS + 3); i++) {
            for (j = 0; j < (OFDM_NC + 2); j++) {
                rx_sym_log[(OFDM_NS + 3)*f+i][j].real = crealf(ofdm->rx_sym[i][j]);
                rx_sym_log[(OFDM_NS + 3)*f+i][j].imag = cimagf(ofdm->rx_sym[i][j]);
            }
        }

        /* note corrected phase (rx no phase) is one big linear array for frame */

        for (i = 0; i < OFDM_ROWSPERFRAME*OFDM_NC; i++) {
            rx_np_log[OFDM_ROWSPERFRAME*OFDM_NC*f + i].real = crealf(ofdm->rx_np[i]);
            rx_np_log[OFDM_ROWSPERFRAME*OFDM_NC*f + i].imag = cimagf(ofdm->rx_np[i]);
        }

        /* note phase/amp ests the same for each col, but check them all anyway */

        for (i = 0; i < OFDM_ROWSPERFRAME; i++) {
            for (j = 0; j < OFDM_NC; j++) {
                phase_est_pilot_log[OFDM_ROWSPERFRAME*f+i][j] = ofdm->aphase_est_pilot_log[OFDM_NC*i+j];
                rx_amp_log[OFDM_ROWSPERFRAME*OFDM_NC*f+OFDM_NC*i+j] = ofdm->rx_amp[OFDM_NC*i+j];
            }
        }

        foff_hz_log[f] = ofdm->foff_est_hz;
        timing_est_log[f] = ofdm->timing_est + 1;     /* offset by 1 to match Octave */
        timing_valid_log[f] = ofdm->timing_valid;     /* offset by 1 to match Octave */
        timing_mx_log[f] = ofdm->timing_mx;           /* offset by 1 to match Octave */
        coarse_foff_est_hz_log[f] = ofdm->coarse_foff_est_hz;
        sample_point_log[f] = ofdm->sample_point + 1; /* offset by 1 to match Octave */

        memcpy(&rx_bits_log[OFDM_BITSPERFRAME*f], rx_bits, sizeof(rx_bits));
    }

    /*---------------------------------------------------------*\
               Dump logs to Octave file for evaluation
                      by tofdm.m Octave script
    \*---------------------------------------------------------*/

    fout = fopen("tofdm_out.txt","wt");
    assert(fout != NULL);
    fprintf(fout, "# Created by tofdm.c\n");
    octave_save_complex(fout, "W_c", (COMP*)ofdm->W, OFDM_NC + 2, OFDM_M, OFDM_M);
    octave_save_complex(fout, "pilot_samples_c", (COMP*)ofdm->pilot_samples, 1, OFDM_M+OFDM_NCP, OFDM_M+OFDM_NCP);
    octave_save_int(fout, "tx_bits_log_c", tx_bits_log, 1, OFDM_BITSPERFRAME*NFRAMES);
    octave_save_complex(fout, "tx_log_c", (COMP*)tx_log, 1, samples_per_frame*NFRAMES,  samples_per_frame*NFRAMES);
    octave_save_complex(fout, "rx_log_c", (COMP*)rx_log, 1, samples_per_frame*NFRAMES,  samples_per_frame*NFRAMES);
    octave_save_complex(fout, "rxbuf_in_log_c", (COMP*)rxbuf_in_log, 1, nin_tot, nin_tot);
    octave_save_complex(fout, "rxbuf_log_c", (COMP*)rxbuf_log, 1, OFDM_RXBUF*NFRAMES,  OFDM_RXBUF*NFRAMES);
    octave_save_complex(fout, "rx_sym_log_c", (COMP*)rx_sym_log, (OFDM_NS + 3)*NFRAMES, OFDM_NC + 2, OFDM_NC + 2);
    octave_save_float(fout, "phase_est_pilot_log_c", (float*)phase_est_pilot_log, OFDM_ROWSPERFRAME*NFRAMES, OFDM_NC, OFDM_NC);
    octave_save_float(fout, "rx_amp_log_c", (float*)rx_amp_log, 1, OFDM_ROWSPERFRAME*OFDM_NC*NFRAMES, OFDM_ROWSPERFRAME*OFDM_NC*NFRAMES);
    octave_save_float(fout, "foff_hz_log_c", foff_hz_log, NFRAMES, 1, 1);
    octave_save_int(fout, "timing_est_log_c", timing_est_log, NFRAMES, 1);
    octave_save_int(fout, "timing_valid_log_c", timing_valid_log, NFRAMES, 1);
    octave_save_float(fout, "timing_mx_log_c", timing_mx_log, NFRAMES, 1, 1);
    octave_save_float(fout, "coarse_foff_est_hz_log_c", coarse_foff_est_hz_log, NFRAMES, 1, 1);
    octave_save_int(fout, "sample_point_log_c", sample_point_log, NFRAMES, 1);
    octave_save_complex(fout, "rx_np_log_c", (COMP*)rx_np_log, 1, OFDM_ROWSPERFRAME*OFDM_NC*NFRAMES, OFDM_ROWSPERFRAME*OFDM_NC*NFRAMES);
    octave_save_int(fout, "rx_bits_log_c", rx_bits_log, 1, OFDM_BITSPERFRAME*NFRAMES);
    fclose(fout);

    ofdm_destroy(ofdm);

    return 0;
}

