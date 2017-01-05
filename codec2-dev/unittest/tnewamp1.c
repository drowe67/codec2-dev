/*---------------------------------------------------------------------------*\

  FILE........: tnewamp1.c
  AUTHOR......: David Rowe
  DATE CREATED: Jan 2017

  Tests for the C version of the newamp1 amplitude modelling used for
  700c.  This program outputs a file of Octave vectors that are loaded
  and automatically tested against the Octave version of the modem by
  the Octave script tnewamp1.m

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

#include "defines.h"
#include "codec2_fft.h"
#include "sine.h"
#include "nlp.h"
#include "dump.h"
#include "octave.h"
#include "quantise.h"

#define FRAMES 100

int main(int argc, char *argv[]) {
    short buf[N_SAMP];	        /* input/output buffer                   */
    float Sn[M_PITCH];	        /* float input speech samples            */
    COMP  Sw[FFT_ENC];	        /* DFT of Sn[]                           */
    codec2_fft_cfg fft_fwd_cfg; /* fwd FFT states                        */
    float w[M_PITCH];	        /* time domain hamming window            */
    COMP  W[FFT_ENC];	        /* DFT of w[]                            */
    MODEL model;
    void *nlp_states;
    float pitch, prev_uq_Wo;
    int   i,m,f;
    
    if (argc != 2) {
        printf("usage: ./tnewamp1 RawFile\n");
        exit(1);
    }
    nlp_states = nlp_create(M_PITCH);
    prev_uq_Wo = TWO_PI/P_MAX;
    fft_fwd_cfg = codec2_fft_alloc(FFT_ENC, 0, NULL, NULL); 
    make_analysis_window(fft_fwd_cfg, w, W);

    for(i=0; i<M_PITCH; i++) {
	Sn[i] = 1.0;
    }

    int K = 20;
    float rate_K_sample_freqs_kHz[K];
    float model_octave[FRAMES][MAX_AMP+2]; // model params in matrix format, useful for C <-> Octave  
    float rate_K_surface[FRAMES][K];       // rate K vecs for each frame, form a surface that makes pretty graphs

    for(f=0; f<FRAMES; f++)
        for(m=0; m<MAX_AMP+2; m++)
            model_octave[f][m] = 0.0;

    mel_sample_freqs_kHz(rate_K_sample_freqs_kHz, K);
    //for(int k=0; k<K; k++)
    //    printf("k: %d sf: %f\n", k, rate_K_sample_freqs_kHz[k]);

    FILE *fin = fopen(argv[1], "rb");
    if (fin == NULL) {
        fprintf(stderr, "Problem opening hts1.raw\n");
        exit(1);
    }

    for(f=0; f<FRAMES; f++) {
        assert(fread(buf,sizeof(short),N_SAMP,fin) == N_SAMP);

        /* shift buffer of input samples, and insert new samples */

	for(i=0; i<M_PITCH-N_SAMP; i++) {
	    Sn[i] = Sn[i+N_SAMP];
	}
	for(i=0; i<N_SAMP; i++) {
	    Sn[i+M_PITCH-N_SAMP] = buf[i];
        }

	/* Estimate Sinusoidal Model Parameters ----------------------*/

	nlp(nlp_states, Sn, N_SAMP, P_MIN, P_MAX, &pitch, Sw, W, &prev_uq_Wo);
	model.Wo = TWO_PI/pitch;

	dft_speech(fft_fwd_cfg, Sw, Sn, w);
	two_stage_pitch_refinement(&model, Sw);
	estimate_amplitudes(&model, Sw, W, 1);

        /* Resample at rate K ----------------------------------------*/

        resample_const_rate_f(&model, &rate_K_surface[f][0], rate_K_sample_freqs_kHz, K);
        for(int k=0; k<K; k++)
            printf("k: %d sf: %f sv: %f\n", k, rate_K_sample_freqs_kHz[k], rate_K_surface[f][k]);
        printf("\n");

        /* log vectors */
 
        model_octave[f][0] = model.Wo;
        model_octave[f][1] = model.L;
        for(m=1; m<=model.L; m++) {
            model_octave[f][m+1] = model.A[m];
        }        
     }

    fclose(fin);

    /* save vectors in Octave format */

    FILE *fout = fopen("tnewamp1_out.txt","wt");
    assert(fout != NULL);
    fprintf(fout, "# Created by tnewamp1.c\n");
    octave_save_float(fout, "rate_K_surface_c", (float*)rate_K_surface, FRAMES, K, K);
    octave_save_float(fout, "model_c", (float*)model_octave, FRAMES, MAX_AMP+2, MAX_AMP+2);
    fclose(fout);

    printf("Done! Now run\n  octave:1> tnewamp1(\"../build_linux/src/hts1a\")\n");
    return 0;
}

