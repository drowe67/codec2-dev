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
#include "newamp1.h"
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
    codec2_fft_cfg phase_fft_fwd_cfg, phase_fft_inv_cfg;
    float pitch, prev_uq_Wo;
    int   i,m,f,k;
    
    if (argc != 2) {
        printf("usage: ./tnewamp1 RawFile\n");
        exit(1);
    }
    nlp_states = nlp_create(M_PITCH);
    prev_uq_Wo = TWO_PI/P_MAX;
    fft_fwd_cfg = codec2_fft_alloc(FFT_ENC, 0, NULL, NULL); 
    make_analysis_window(fft_fwd_cfg, w, W);

    phase_fft_fwd_cfg = codec2_fft_alloc(NEWAMP1_PHASE_NFFT, 0, NULL, NULL);
    phase_fft_inv_cfg = codec2_fft_alloc(NEWAMP1_PHASE_NFFT, 1, NULL, NULL);

    for(i=0; i<M_PITCH; i++) {
	Sn[i] = 1.0;
    }

    int K = 20;
    float rate_K_sample_freqs_kHz[K];
    float model_octave[FRAMES][MAX_AMP+2];    // model params in matrix format, useful for C <-> Octave  
    float rate_K_surface[FRAMES][K];          // rate K vecs for each frame, form a surface that makes pretty graphs
    float rate_K_surface_no_mean[FRAMES][K];  // mean removed surface  
    float rate_K_surface_no_mean_[FRAMES][K]; // quantised mean removed surface  
    float mean[FRAMES];
    float mean_[FRAMES];
    float rate_K_surface_[FRAMES][K];         // quantised rate K vecs for each frame
    float interpolated_surface_[FRAMES][K];   // dec/interpolated surface
    int   voicing[FRAMES];
    int   voicing_[FRAMES];
    float model_octave_[FRAMES][MAX_AMP+2];
    COMP  H[FRAMES][MAX_AMP];
    int indexes[FRAMES][NEWAMP1_N_INDEXES];

    for(f=0; f<FRAMES; f++) {
        for(m=0; m<MAX_AMP+2; m++) {
            model_octave[f][m] = 0.0;
            model_octave_[f][m] = 0.0;
        }
        for(m=0; m<MAX_AMP; m++) {
            H[f][m].real = 0.0;
            H[f][m].imag = 0.0;
        }
        for(k=0; m<K; k++)
            interpolated_surface_[f][k] = 0.0;
        voicing_[f] = 0;
    }

    mel_sample_freqs_kHz(rate_K_sample_freqs_kHz, K);

    //for(int k=0; k<K; k++)
    //    printf("k: %d sf: %f\n", k, rate_K_sample_freqs_kHz[k]);

    FILE *fin = fopen(argv[1], "rb");
    if (fin == NULL) {
        fprintf(stderr, "Problem opening hts1.raw\n");
        exit(1);
    }

    int M = 4; 

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
        est_voicing_mbe(&model, Sw, W);
        voicing[f] = model.voiced;

        /* newamp1 processing ----------------------------------------*/

        newamp1_model_to_indexes(&indexes[f][0], 
                                 &model, 
                                 &rate_K_surface[f][0], 
                                 rate_K_sample_freqs_kHz,
                                 K,
                                 &mean[f],
                                 &rate_K_surface_no_mean[f][0],
                                 &rate_K_surface_no_mean_[f][0]);

        newamp1_indexes_to_rate_K_vec(&rate_K_surface_[f][0],
                                      &rate_K_surface_no_mean_[f][0],
                                      rate_K_sample_freqs_kHz,
                                      K,
                                      &mean_[f],
                                      &indexes[f][0]);

        /* log vectors */
 
        model_octave[f][0] = model.Wo;
        model_octave[f][1] = model.L;
        for(m=1; m<=model.L; m++) {
            model_octave[f][m+1] = model.A[m];
        }        
    }

 
    /* Decoder */

    MODEL model__[M];
    float prev_rate_K_vec_[K];
    COMP  HH[M][MAX_AMP+1];
    float Wo_left;
    int   voicing_left;

    for(k=0; k<K; k++)
        prev_rate_K_vec_[k] = rate_K_surface_[0][k];

    for(f=0; f<FRAMES; f+=M) {

#ifdef TMP1
        /* Quantise Wo. V/UV flag is coded using a zero index for Wo,
           this means we need to adjust Wo index slightly for the
           lowest Wo V frames */

        if (voicing[f]) {
            int index = encode_log_Wo(model_octave[f][0], 6);
            if (index == 0) {
                index = 1;
            }
            model_octave_[f][0] = decode_log_Wo(index, 6);
         }
        else {
            model_octave_[f][0] = 2.0*M_PI/100.0;
        }
#endif

        if (indexes[f][3]) {
            model_octave_[f][0] = decode_log_Wo(indexes[f][3], 6);
        }
        else {
            model_octave_[f][0] = 2.0*M_PI/100.0;
        }

        /*
          pass in left and right rate K vectors, Wo, v
          interpolate at rate K
          convert to rate L
          indexes for current plus decoded rate K, Wo, v from frame[0], model[4] out parameters,
          then plug them into model_octave for testing
          ref indexes to 0...3
          slowly change
          
          [ ] model[4] out
          [ ] change to 1:M processing
              + not sure how to handle Octave side of this
        */

        if (f >= M) {
            /*
            float *left_vec = &rate_K_surface_[f-M][0];
            float *right_vec = &rate_K_surface_[f][0];
            newamp1_interpolate(&interpolated_surface_[f-M][0], left_vec, right_vec, K);
            */

            Wo_left = model_octave_[f-M][0];
            voicing_left = voicing[f-M];

            newamp1_indexes_to_model(model__,
                                     (COMP*)HH,
                                     &interpolated_surface_[f-M][0],
                                     prev_rate_K_vec_,
                                     &Wo_left,
                                     &voicing_left,
                                     rate_K_sample_freqs_kHz, 
                                     K,
                                     phase_fft_fwd_cfg, 
                                     phase_fft_inv_cfg,
                                     &indexes[f][0]);

            /* interpolate 25Hz v and Wo back to 100Hz */

            float aWo_[M];
            int avoicing_[M], aL_[M], m;
            float Wo1 = model_octave_[f-M][0];
            float Wo2 = model_octave_[f][0];
            interp_Wo_v(aWo_, aL_, avoicing_, Wo1, Wo2, voicing[f-M], voicing[f]);

            for(i=f-M, m=0; i<f; i++,m++) {
                model_octave_[i][0] = model__[m].Wo; //aWo_[m];
                model_octave_[i][1] = model__[m].L; //aL_[m]; 
                voicing_[i] = avoicing_[m];
            }

            /* back to rate L, synth phase */

            MODEL model_;
            COMP aH[MAX_AMP+1];
            int j;
            for(i=f-M, j=0; i<f; i++,j++) {
                /*
                model_.Wo = model_octave_[i][0];
                model_.L  = model_octave_[i][1];
                resample_rate_L(&model_, &interpolated_surface_[i][0], rate_K_sample_freqs_kHz, K);
                determine_phase(aH, &model_, NEWAMP1_PHASE_NFFT, phase_fft_fwd_cfg, phase_fft_inv_cfg);
                */
                model_.L  = model_octave_[i][1];
                for(m=1; m<=model_.L; m++) {
                    model_octave_[i][m+1] = model__[j].A[m]; // = model_.A[m];
                    H[i][m-1] = HH[j][m];// aH[m];
                }
            }
        }
    }

    fclose(fin);

    /* save vectors in Octave format */

    FILE *fout = fopen("tnewamp1_out.txt","wt");
    assert(fout != NULL);
    fprintf(fout, "# Created by tnewamp1.c\n");
    octave_save_float(fout, "rate_K_surface_c", (float*)rate_K_surface, FRAMES, K, K);
    octave_save_float(fout, "mean_c", (float*)mean, 1, FRAMES, 1);
    octave_save_float(fout, "rate_K_surface_no_mean_c", (float*)rate_K_surface_no_mean, FRAMES, K, K);
    octave_save_float(fout, "rate_K_surface_no_mean__c", (float*)rate_K_surface_no_mean_, FRAMES, K, K);
    octave_save_float(fout, "mean__c", (float*)mean_, FRAMES, 1, 1);
    octave_save_float(fout, "rate_K_surface__c", (float*)rate_K_surface_, FRAMES, K, K);
    octave_save_float(fout, "interpolated_surface__c", (float*)interpolated_surface_, FRAMES, K, K);
    octave_save_float(fout, "model_c", (float*)model_octave, FRAMES, MAX_AMP+2, MAX_AMP+2);
    octave_save_float(fout, "model__c", (float*)model_octave_, FRAMES, MAX_AMP+2, MAX_AMP+2);
    octave_save_int(fout, "voicing__c", (int*)voicing_, 1, FRAMES);
    octave_save_complex(fout, "H_c", (COMP*)H, FRAMES, MAX_AMP, MAX_AMP);
    fclose(fout);

    printf("Done! Now run\n  octave:1> tnewamp1(\"../build_linux/src/hts1a\")\n");
    return 0;
}

