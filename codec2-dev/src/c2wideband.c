/*---------------------------------------------------------------------------*\

  FILE........: c2wideband.c
  AUTHOR......: David Rowe & Phil Ayres
  DATE CREATED: July 2017

 * TODO - what is this file doing?
 * 
\*---------------------------------------------------------------------------*/

/*
  Copyright David Rowe 2017

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

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>

#include "defines.h"
#include "codec2_fft.h"
#include "sine.h"
#include "nlp.h"
#include "dump.h"
#include "lpc.h"
#include "quantise.h"
#include "phase.h"
#include "interp.h"
#include "postfilter.h"
#include "codec2.h"
#include "lsp.h"
#include "codec2_internal.h"
#include "machdep.h"
#include "codec2.h"
#include "newamp1.h"
#include "dct2.h"

#include "c2wideband.h"
#include "c2wideband_map.h"

float mean(float data[], int n);
int unit_test();
float std(float data[], int rows);
void diff_de(int rows, int cols, float D[rows][cols], float E[rows][cols], float diff_de[rows][cols]);
void array_col_to_row(int rows, int cols, float data[rows][cols], int col, float res[]);
void std_on_cols(int rows, int cols, float data[rows][cols], float res[]);
void setup_map(WIDEBAND_MAP * wb_map, int Nt, int K);

#define C2WB_PLOT


//TODO: move all this to the standard codec2.c functions
#ifdef IGNORE

int main(int argc, char *argv[])
{
    struct CODEC2 *codec2;
    FILE *fin;
    FILE *fout;
    short *buf;
    unsigned char *bits;
    int nsam, nbit, i, r;
    int errno = -1;
    for (i = 0; i < 10; i++) {
        r = codec2_rand();
        printf("[%d] r = %d\n", i, r);
    }

    if (argc != 3) {
        printf("usage: %s InputRawSpeechFile OutputRawSpeechFile\n", argv[0]);
        exit(1);
    }

    if ((fin = fopen(argv[1], "rb")) == NULL) {
        fprintf(stderr, "Error opening input speech file: %s: %s.\n",
                argv[1], strerror(errno));
        exit(1);
    }

    if ((fout = fopen(argv[2], "wb")) == NULL) {
        fprintf(stderr, "Error opening output speech file: %s: %s.\n",
                argv[2], strerror(errno));
        exit(1);
    }

#ifdef DUMP
    dump_on("c2demo");
#endif

    /* Note only one set of Codec 2 states is required for an encoder
       and decoder pair. */

    codec2 = codec2_create(CODEC2_MODE_WB);
    nsam = C2WB_SPERF; //codec2_samples_per_frame(codec2);
    buf = (short*) malloc(nsam * sizeof (short));
    nbit = C2WB_BPERF; //codec2_bits_per_frame(codec2);
    bits = (unsigned char*) malloc(nbit * sizeof (char));

    while (fread(buf, sizeof (short), nsam, fin) == (size_t) nsam) {
        codec2_encode_wb(codec2, bits, buf);
        codec2_decode_wb(codec2, buf, bits);
        fwrite(buf, sizeof (short), nsam, fout);
    }

    free(buf);
    free(bits);
    codec2_destroy(codec2);

    fclose(fin);
    fclose(fout);

    return 0;
}
#endif

void codec2_encode_wb(struct CODEC2 *c2, unsigned char * bits, short speech[])
{
    return;
}

void codec2_decode_wb(struct CODEC2 *c2, short speech[], const unsigned char * bits)
{
    return;
}

void calculate_Am_freqs_kHz(float Wo, int L, float Am_freqs_kHz[L])
{
    //(1:L)*Wo*4/pi;
    int i;
    for (i = 0; i < L; i++) {
        Am_freqs_kHz[i] = i * Wo * 4 / PI;
    }
}


// Adapted from newamp.m

void resample_const_rate_f_mel(C2CONST *c2const, MODEL * model, float K, float* rate_K_surface, float* rate_K_sample_freqs_kHz)
{
    mel_sample_freqs_kHz(rate_K_sample_freqs_kHz, K, 100, 0.95 * c2const->Fs / 2);
    resample_const_rate_f(c2const, model, rate_K_surface, rate_K_sample_freqs_kHz, K);
}


// Updates rate_K_vec with the corrected values
// function [rate_K_vec_corrected orig_error error nasty_error_log nasty_error_m_log] = correct_rate_K_vec(rate_K_vec, rate_K_sample_freqs_kHz, AmdB, AmdB_, K, Wo, L, Fs)

void correct_rate_K_vec(MODEL *model, float rate_K_vec[], float rate_K_sample_freqs_kHz[], float Am_freqs_kHz[], float orig_AmdB[], int K, float Wo, int L, int Fs, float rate_K_vec_corrected[])
{

    /*
        % aliasing correction --------------------------------------

        % The mel sample rate decreases as frequency increases. Look for
        % any regions above 1000Hz where we have missed definition of a
        % spectral peak (formant) due to aliasing.  Adjust the rate K
        % sample levels to restore peaks.  Theory is that correct
        % definition of a formant is less important than the frequency of
        % the formant.  As long as we define a formant in that general
        % frequency area it will sound OK.
     */


    // this is passed in
    //Am_freqs_kHz = (1:L)*Wo*4/pi;    

    //% Lets see where we have made an error
    float error[MAX_AMP + 1];
    //float orig_error[MAX_AMP+1];
    float AmdB_[MAX_AMP + 1];

    float nasty_error_freq;
    float rate_K_prev_sample_kHz;
    float rate_K_next_sample_kHz;

    int Ncorrections = 3; //% maximum number of rate K samples to correct
    int error_thresh = 3; //% only worry about errors larger than thresh    

    int closest_k = 0;
    int i;

    // regenerate the AmdB values from the updated model
    for (i = 0; i < MAX_AMP + 1; i++) {
        AmdB_[i] = 20 * log10(model->A[i]);
    }

    // calculate error between original AmdB and the new values
    for (i = 0; i < MAX_AMP + 1; i++) {
        int a = orig_AmdB[i] - AmdB_[i];
        error[i] = a;
        //  orig_error[i] = a;
    }



    //% first 1000Hz is densely sampled so ignore
    int start_m = floor(L * 1000 / (Fs / 2));
    for (i = 0; i < start_m; i++) {
        error[i] = 0;
    }
    //start_m = floor(L*1000/(Fs/2));
    //error(1:start_m) = 0;  % first 1000Hz is densly sampled so ignore

    // ignore these
    //float nasty_error_m_log = []; 
    //float nasty_error_log = [];


    // could probably memcpy this, but I need to think about it...
    // TODO copy??
    for (i = 0; i < K; i++) {
        //rate_K_vec_corrected = rate_K_vec
        rate_K_vec_corrected[i] = rate_K_vec[i];
    }

    for (i = 0; i < Ncorrections; i++) {

        // [mx mx_m] = max(error);
        int mx_m = 0;
        int mx = 0;
        int m;
        for (m = 0; m < start_m; m++) {
            if (error[m] > mx) {
                mx_m = m;
                mx = error[m];
            }
        }

        if (mx > error_thresh) {

            // ignore these
            //nasty_error_log = [nasty_error_log mx];
            //nasty_error_m_log = [nasty_error_m_log mx_m];

            //% find closest rate K sample to nasty error

            nasty_error_freq = (float) mx_m * Wo * Fs / (2 * PI * 1000);

            //[tmp closest_k] = min(abs(rate_K_sample_freqs_kHz - nasty_error_freq));
            closest_k = -1;
            float cka = K;
            int k;
            for (k = 0; k < K; k++) {
                float a = abs(rate_K_sample_freqs_kHz[k] - nasty_error_freq);

                if (closest_k == -1 || a < cka) {
                    closest_k = k;
                    cka = a;
                }
            }

            // this sets the value from the original AmdB value set for the model prior to recalculation
            // TODO - check this is correct
            rate_K_vec_corrected[closest_k] = orig_AmdB[mx_m];

            //% zero out error in this region and look for another large error region

            // here I'm assuming that the required result is a zero based index
            // rather than an Octave 1 based index
            // closest_k is already 0 based in this code
            // TODO - check this
            k = (0 > closest_k - 1 ? 0 : closest_k - 1);
            //k = max(1, closest_k-1); 

            rate_K_prev_sample_kHz = rate_K_sample_freqs_kHz[k];


            // again, assuming the required result is 0 rather than 1 based
            // TODO - check 
            k = (K - 1 < closest_k + 1 ? K : closest_k + 1);
            //k = min(K, closest_k+1); 

            rate_K_next_sample_kHz = rate_K_sample_freqs_kHz[k];

            int st = -1;
            int st_m = 0;
            int en = -1;
            int en_m = 0;
            for (m = 0; m < C2WB_K - 1; m++) {
                //[tmp st_m] = min(abs(Am_freqs_kHz - rate_K_prev_sample_kHz));
                //[tmp en_m] = min(abs(Am_freqs_kHz - rate_K_next_sample_kHz));
                int pa = abs(Am_freqs_kHz[m] - rate_K_prev_sample_kHz);
                int na = abs(abs(Am_freqs_kHz[m] - rate_K_next_sample_kHz));
                if (st == -1 || pa < st) {
                    st = pa;
                    st_m = m;
                }
                if (en == -1 || na < en) {
                    en = na;
                    en_m = m;
                }
            }
            if (closest_k == K)
                en_m = L;

            for (m = st_m; m < en_m; m++) {
                error[m] = 0;
            }
            //error(st_m:en_m) = 0;
        }
    }



}

void diff_de(int rows, int cols, float D[rows][cols], float E[rows][cols], float diff_de[rows][cols])
{

    int col, row;
    for (col = 0; col < cols; col++) {
        for (row = 0; row < rows; row++) {
            diff_de[row][col] = D[row][col] - E[row][col];
        }
    }

}

float mean(float data[], int n)
{
    float sum = 0.0;
    int i;
    for (i = 0; i < n; ++i) {
        sum += data[i];
    }

    return sum / n;
}

void array_col_to_row(int rows, int cols, float data[rows][cols], int col, float res[])
{
    int row;
    for (row = 0; row < rows; row++) {
        res[row] = data[row][col];
    }
}

float std(float data[], int rows)
{
    float standardDeviation = 0.0;
    int i;
    for (i = 0; i < rows; ++i) {
        standardDeviation += pow(data[i] - mean(data, rows), 2);
    }
    return sqrt(standardDeviation / rows);
}

void std_on_cols(int rows, int cols, float data[rows][cols], float res[])
{
    float row_from_col[cols];
    int col;
    for (col = 0; col < cols; col++) {
        array_col_to_row(rows, cols, data, col, row_from_col);
        res[col] = std(row_from_col, rows);
    }
}

float mean_std_diff_de(int rows, int cols, float D[rows][cols], float E[rows][cols])
{
    float matrix_diff_de[rows][cols];
    float std_array[cols];
    diff_de(rows, cols, D, E, matrix_diff_de);
    std_on_cols(rows, cols, matrix_diff_de, std_array);
    return mean(std_array, cols);
}


//% Encode/decoder a 160ms block of model parameters
//% TODO: (i) quantisation of DCT coeffs (ii) break into separate encoder and decoder functions
//
//function [model_block_ dct2_sd qn rate_K_surface_block rate_K_surface_block_] = wideband_enc_dec(model_block, rmap, cmap)

void wideband_enc_dec(C2CONST *c2const, int n_block_frames, MODEL model_block[n_block_frames], WIDEBAND_MAP * wb_map,
                      MODEL model_block_[n_block_frames], float dct2_sd[n_block_frames], int * p_qn, float rate_K_surface_block[n_block_frames][C2WB_K], float rate_K_surface_block_[n_block_frames][C2WB_K])
{


    //printf("starting wideband_enc_dec with model_block[0]->L = %d\n", model_block[0].L);
    //    c2wideband_const;
    //
    //    sim_quant = 1; % used to simulate quantisation, set to 1,2,4, etc
    //    dist_dB   = 2; % use enough coefficients to get this distortion ond DCT coeffs
    int sim_quant = 1;
    float dist_dB = 2;

    //    int Fs = c2const->Fs; 
    //    int L;   
    //    float Wo;

    int K = C2WB_K;
    int Nt = C2WB_NT;

    int rows = Nt;
    int cols = K;
    int dec = C2WB_DEC;

    // one time configuration of DCT & IDCT

    codec2_dct_cfg dct_cfg_n = dct_config(cols);
    codec2_dct_cfg dct_cfg_m = dct_config(rows);
    codec2_dct_cfg idct_cfg_n = idct_config(cols);
    codec2_dct_cfg idct_cfg_m = idct_config(rows);

    //printf("starting iteration\n");
    // iterate through the frames in the block
    int f;
    for (f = 0; f < n_block_frames; f++) {
        //printf("processing frame %d of %d\n", f, n_block_frames);

        float rate_K_sample_freqs_kHz[C2WB_K];
        //MODEL * model = &model_block[f];


        //    % Resample variable rate L vectors to fixed length rate K.  We have
        //    % left high end correction out for now, this is less of an issue
        //    % with a higher K


        //    [rate_K_surface_block rate_K_sample_freqs_kHz] = resample_const_rate_f_mel(model_block, K, Fs);
        //printf("resample_const_rate_f_mel\n");
        resample_const_rate_f_mel(c2const, &model_block[f], K, rate_K_surface_block[f], rate_K_sample_freqs_kHz);

        //    % decimate down to 20ms time resolution, and DCT

        //    D = dct2(rate_K_surface_block(1:dec:Nt*dec,:));

        float D[rows][cols];
        float E[rows][cols];
        //printf("dct2\n");
        dct2(dct_cfg_m, dct_cfg_n, rows, cols, &rate_K_surface_block[f], D);

        //    % So D is the 2D block of DCT coeffs at the encoder.  We want to
        //    % create a quantised version at the "decoder" E.  This loop copies
        //    % DCTs coeffs from D to E, until we get beneath a distortion
        //    % threshold.
        //
        //    % This is essentially variable rate quantisation, but gives us
        //    % some idea of the final bit rate.  In practice we will also need
        //    % to constrain the total number of bits (ie bit rate), and
        //    % quantise each coefficient.
        //
        //    % Turns out than mean SD (across many blocks/frames) is about the
        //    % same in the DCT domain as the rate K domain.  So we can just
        //    % measure MSE between D and E to estimate mean SD of the rate K
        //    % vectors after quantisation.

        //    E = mapped = zeros(Nt,K);                
        int r, c;
        for (r = 0; r < rows; r++) {
            memset(E[r], '\0', cols * sizeof (float));
        }

        //    qn = 0;
        *p_qn = 0;

        //    adct2_sd = mean(std(D-E));

        float adct2_sd;


        adct2_sd = mean_std_diff_de(rows, cols, D, E);

        //printf("while adct2_sd > dist_dB (dist_dB=%f)\n", dist_dB);
        //    while adct2_sd > dist_dB
        //TODO : validate the addition of  *p_qn limit check is correct
        while (adct2_sd > dist_dB && *p_qn < C2WB_NT) {

            //      qn++;                        
            // --- this has moved to the end to cope with 0 base
            assert(*p_qn < C2WB_NT);
            int rmapqn = wb_map->rmap[*p_qn];
            int cmapqn = wb_map->cmap[*p_qn];
            assert(rmapqn < rows);
            assert(cmapqn < cols);
            //      E(rmap(qn), cmap(qn)) = sim_quant*round(D(rmap(qn), cmap(qn))/sim_quant);
            E[rmapqn][cmapqn] = sim_quant * round(D[rmapqn][cmapqn] / sim_quant);
            //      adct2_sd = mean(std(D-E));
            adct2_sd = mean_std_diff_de(rows, cols, D, E);

            //printf("qn %d %f\n", *p_qn, adct2_sd);

            (*p_qn)++;
            //      %printf("qn %d %f\n", qn, adct2_sd);
            //    end
        }

        //    % note neat trick to interpolate to 10ms frames despite dec to 20ms, this means
        //    % we don't need a separate decode side interpolator.

        //    dct2_sd = mean(std(D-E));


        dct2_sd[f] = mean_std_diff_de(rows, cols, D, E);


        //    rate_K_surface_block_ = idct2([sqrt(dec)*E; zeros(Nt*(dec-1), K)]);
        float inrks[rows * dec][K];

        //printf("setting inrks\n");
        for (r = 0; r < rows; r++) {
            for (c = 0; c < cols; c++) {
                inrks[r][c] = sqrt(dec) * E[r][c];
            }
        }
        for (r = 0; r < Nt * (dec - 1); r++) {
            for (c = 0; c < K; c++) {
                inrks[r + rows][c] = 0;
            }
        }


        //TODO ???
        //= [sqrt(dec)*E; zeros(Nt*(dec-1), K)];
        //printf("idct2\n");
        idct2(idct_cfg_m, idct_cfg_n, rows, cols, inrks, &rate_K_surface_block_[f]);

        //    model_block_ = resample_rate_L(model_block, rate_K_surface_block_, rate_K_sample_freqs_kHz, Fs);        
        //printf("resample_rate_L\n");
        resample_rate_L(c2const, &model_block[f], rate_K_surface_block_[f], rate_K_sample_freqs_kHz, K);
        //endfunction 
    }
}

void setup_map(WIDEBAND_MAP * wb_map, int Nt, int K)
{

    /*
      % map that defines order we read out and quantise DCT coeffs
      % TODO: for C port we need an Octave function to write Map to a C
      % include file

      map = load("c2wideband_map");

      % create arrays to reverse map quantiser_num to r,c Luts

      rmap = cmap = zeros(1,Nt*K);
      for r=1:Nt
        for c=1:K
           quantiser_num = map(r,c);
           rmap(quantiser_num) = r;
           cmap(quantiser_num) = c;
        end
      end

     */

    //printf("setup_map(wb_map, Nt, K)");
    int quantiser_num, r, c;

    memset(wb_map->rmap, '\0', Nt * K * sizeof *wb_map->rmap);
    memset(wb_map->cmap, '\0', Nt * K * sizeof *wb_map->cmap);

    for (r = 0; r < Nt; r++) {
        for (c = 0; c < K; c++) {
            quantiser_num = c2wideband_map[r][c];
            wb_map->rmap[quantiser_num] = r;
            wb_map->cmap[quantiser_num] = c;
        }
    }

}

/*
% ---------------------------------------------------------------------------------------
% rate K mel-resampling, high end correction, and DCT experiment workhorse

function [model_ rate_K_surface] = experiment_rate_K_dct2(model, plots=1)
 */

//TODO - review and produce a unit test for this
// this is really just a draft at this point

void experiment_rate_K_dct2(C2CONST *c2const, MODEL model_frames[], const int total_frames)
{

    printf("experiment_rate_K_dct2 with frames: %d\n", total_frames);
    //  newamp;
    //  c2wideband_const;
    //  [frames nc] = size(model);

    const int n_block_frames = C2WB_NT * C2WB_DEC;

    const int K = C2WB_K;
    const int Nt = C2WB_NT;
    const int dec = C2WB_DEC;
    const int Tf = C2WB_TF;
    WIDEBAND_MAP wb_map;


    //  % break into blocks of (Nt time samples) x (K freq samples)

    // Nblocks = floor(frames/(Nt*dec));

    const int Nblocks = floor(total_frames / n_block_frames);

    // number of frames that can be processed (if the final block is not a full set of frames)
    int frames;

    frames = Nblocks * n_block_frames;

    printf("total_frames: %d processable frames: %d Nblocks: %d\n", total_frames, frames, Nblocks);


    setup_map(&wb_map, Nt, K);



    //% per-block processing ----------------------------------------------------

    //  % init a bunch of output variables

    //  rate_K_surface_ = zeros(Nblocks*Nt*dec, K);  
    // sumnz = zeros(1,Nblocks);
    // dct2_sd = zeros(1,Nblocks);
    // model_ = [];


    float rate_K_surface_block[total_frames][K]; // rate K vecs for each frame, form a surface that makes pretty graphs   
    float rate_K_surface_block_[total_frames][K];

#ifdef C2WB_PLOT    
    //  float rate_K_surface[total_frames][K];
    //  float rate_K_surface_[total_frames][K];
    //  MODEL model_[total_frames];
#endif

    float sumnz[Nblocks];
    float dct2_sd[Nblocks];
    int qn;

    MODEL model_block_[n_block_frames];

    //for n=1:Nblocks
    // Step through the model in blocks of 16 frames
    int n_block = 0;
    int f;
    for (f = 0; f < frames; f += n_block_frames) {
        MODEL * model_block = &model_frames[f];


        //st = (n-1)*dec*Nt+1; en = st + dec*Nt - 1;
        //printf("st: %d en: %d\n", st, en);
        // effectively handled by the iterations through the block in wideband_enc_dec

        //[model_block_ adct2_sd qn rate_K_surface_block rate_K_surface_block_] = wideband_enc_dec(model(st:en,:), rmap, cmap);
        //    void wideband_enc_dec(C2CONST *c2const, int n_block_frames, MODEL model_block[n_block_frames], WIDEBAND_MAP * wb_map,
        //        MODEL model_block_[n_block_frames], float dct2_sd[n_block_frames], int * p_qn, float rate_K_surface_block[n_block_frames][C2WB_K], float rate_K_surface_block_[n_block_frames][C2WB_K]) {

        wideband_enc_dec(c2const, n_block_frames, model_block, &wb_map,
                         model_block_, &dct2_sd[n_block], &qn, &rate_K_surface_block[f], &rate_K_surface_block_[f]);

#ifdef C2WB_PLOT
        //model_ = [model_; model_block_];
        //% log these for plotting/development

        //rate_K_surface(st:en,:) = rate_K_surface_block;
        //rate_K_surface_(st:en,:) = rate_K_surface_block_;

        //    for (int p = 0; p < n_block_frames; p++) {

        //      model_[f+p] = model_block_[p];

        //      for (int k = 0; k < K; k++) {
        //        rate_K_surface[f+p][k] = rate_K_surface_block[n_block + p][k];
        //        rate_K_surface_[f+p][k] = rate_K_surface_block_[n_block + p][k];
        //      }
        //    }
        //dct2_sd(n) = adct2_sd;            
        sumnz[n_block] = (float) qn;
        // already handled in call
#endif
        n_block++;
        //end
    }


#ifdef C2WB_PLOT
    printf("average dct spectral distortion: %3.2f dB\n", mean(dct2_sd, Nblocks));
    printf("mean number of coeffs/DCT: %3.2f/%d\n", mean(sumnz, Nblocks), Nt * K);
    printf("coeffs/second: %3.2f\n", mean(sumnz, Nblocks) / (Nt * Tf * dec));
    printf("bits/s: %3.2f\n", 2.9 * mean(sumnz, Nblocks) / (Nt * Tf * dec));

    //TODO
    //dist = std((rate_K_surface_(1:dec:Nblocks*Nt*dec,:) - rate_K_surface(1:dec:Nblocks*Nt*dec,:))');
    //float dist = 
    //if plots
    //figure(1); clf; plot(dist); title('Rate K SD');
    //printf("Rate K spectral distortion mean: %3.2f dB var: %3.2f\n", mean(dist), var(dist));
    //end
#endif



}
