#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>

#include "t_helpers.h"
#include "c2wideband.h"
#include "codec2.h"
#include "defines.h"
#include "codec2_fft.h"
#include "sine.h"
#include "nlp.h"
#include "dump.h"
#include "octave.h"
#include "newamp1.h"
#include "quantise.h"

#define FRAMES 160

float mean(float data[], int n);
int unit_test();
float std(float data[], int rows);
void diff_de(int rows, int cols, float D[rows][cols], float E[rows][cols], float diff_de[rows][cols]);
void array_col_to_row(int rows, int cols, float data[rows][cols], int col, float res[]);
void std_on_cols(int rows, int cols, float data[rows][cols], float res[]);
float mean_std_diff_de(int rows, int cols, float D[rows][cols], float E[rows][cols]);
void test_wideband_enc_dec();
void setup_map(WIDEBAND_MAP * wb_map, int Nt, int K);
void test_with_real_data(int argc, char *argv[]);

int main(int argc, char *argv[])
{

    printf("Testing file c2wideband\n");

    test("mean");
    int n = 5;
    float in[] = {1.0, 2.0, 4.0, 6.0, 10.0};
    float res_f;
    float expect_f;

    res_f = mean(in, n);
    expect_f = 23.0 / 5;
    if (res_f != expect_f) {
        test_failed_f(expect_f, res_f);
    }

    test("std");
    res_f = std(in, n);
    expect_f = 3.2000000000000002;
    if (res_f != expect_f) {
        test_failed_f(expect_f, res_f);
    }


    test("diff_de");
    int rows = 5;
    int cols = 3;
    float D[rows][cols];
    float E[rows][cols];
    float res_diff_de[rows][cols];
    float expect_diff_de[rows][cols];
    int r, c;

    for (r = 0; r < rows; r++) {
        for (c = 0; c < cols; c++) {
            float d = rand();
            float e = rand();
            D[r][c] = d;
            E[r][c] = e;
            expect_diff_de[r][c] = d - e;
        }
    }

    diff_de(rows, cols, D, E, res_diff_de);

    if (memcmp(res_diff_de, expect_diff_de, rows * cols * sizeof (float)) != 0) {
        test_failed();
    }


    test("array_col_to_row");
    float data[rows][cols];
    float res_data[rows];
    float expect_data[cols][rows];
    for (r = 0; r < rows; r++) {
        for (c = 0; c < cols; c++) {
            float d = rand();
            data[r][c] = d;
            expect_data[c][r] = d;
        }
    }
    for (c = 0; c < cols; c++) {

        array_col_to_row(rows, cols, data, c, res_data);
        if (memcmp(res_data, expect_data[c], cols * sizeof (float)) != 0) {
            test_failed();
        }
    }


    test("std_on_cols");

    cols = 4;

    float data_std[5][4] = {
        {1.0, 1.0, 4.0, 6.0},
        {1.0, 2.0, 5.5, 6.0},
        {1.0, 4.0, 6.0, 6.0},
        {1.0, 6.0, 3.3, 6.0},
        {1.0, 10.0, -0.2, 6.0}
    };

    float res_data_std[cols];
    float expect_data_std[] = {0.0, 3.2000000000000002, 2.1903424389807178, 0.0};

    std_on_cols(rows, cols, data_std, res_data_std);
    if (memcmp(res_data_std, expect_data_std, cols * sizeof (float)) != 0) {
        test_failed(fn);
    }

    test("mean_std_diff_de");

    cols = 4;

    float data_std_d[5][4] = {
        {2.0, 1.0, 5.0, 6.0},
        {1.0, 2.0, 10.5, 6.0},
        {4.0, 4.0, 12.0, 6.0},
        {1.0, 6.0, 4.6, 6.0},
        {6.5, 10.0, 0.0, 6.0}
    };

    float data_std_e[5][4] = {
        {1.0, 0.0, 1.0, 0.0},
        {0.0, 0.0, 5.0, 0.0},
        {3.0, 0.0, 6.0, 0.0},
        {0.0, 0.0, 1.3, 0.0},
        {5.5, 0.0, 0.2, 0.0}
    };

    float expected_msd = (3.2000000000000002 + 2.1903424389807178) / cols;

    float res_msd;
    res_msd = mean_std_diff_de(rows, cols, data_std_d, data_std_e);
    if (abs(res_msd - expected_msd) > 0.00001) {
        test_failed_f(res_msd, expected_msd);
    }

    test_wideband_enc_dec();

    test_with_real_data(argc, argv);

    return 1;
}

void test_wideband_enc_dec()
{

    test("wideband_enc_dec");

    short *buf;
    unsigned char *bits;
    //int            nsam, nbit, i, r;
    int nbit;
    int nsam;
    int n_block_frames = C2WB_NT * C2WB_DEC;
    int K = C2WB_K;
    int Nt = C2WB_NT;
    int Fs = C2WB_FS;
    int dec = C2WB_DEC;
    int Nblocks = n_block_frames;
    int rate_K_surface_size = Nblocks * Nt * dec;

    struct CODEC2 * codec2 = codec2_create(CODEC2_MODE_WB);

    C2CONST c2const = c2const_create(Fs);
    nsam = C2WB_SPERF; //codec2_samples_per_frame(codec2);
    nbit = C2WB_BPERF; //codec2_bits_per_frame(codec2);

    buf = (short*) malloc(nsam * sizeof (short));
    bits = (unsigned char*) malloc(nbit * sizeof (char));

    WIDEBAND_MAP wb_map;
    setup_map(&wb_map, Nt, K);

    MODEL model_block[n_block_frames];
    MODEL model_block_[n_block_frames];
    float dct2_sd[Nblocks];
    int qn;
    float rate_K_surface_block[rate_K_surface_size][K]; // rate K vecs for each frame, form a surface that makes pretty graphs   
    float rate_K_surface_block_[rate_K_surface_size][K];
    int n, i;

    for (n = 0; n < n_block_frames; n++) {
        model_block[n].L = 10;
        model_block[n].Wo = (float) rand() / RAND_MAX;
        for (i = 0; i < MAX_AMP + 1; i++) {
            model_block[n].phi[i] = (float) rand() / RAND_MAX;
            model_block[n].A[i] = (float) rand() / RAND_MAX;
        }
        model_block[n].voiced = round((float) rand() / RAND_MAX);
    }

    n = 0;

    printf("setup complete. now calling the function\n");

    wideband_enc_dec(&c2const, n_block_frames, model_block, &wb_map,
                     model_block_, dct2_sd, &qn, &rate_K_surface_block[n], &rate_K_surface_block_[n]);

    codec2_destroy(codec2);
    free(bits);
    free(buf);
    printf("made it to the end\n");
}

/*---------------------------------------------------------------------------*\

  FUNCTION....: test_with_real_data()
  AUTHOR......: David Rowe
  DATE CREATED: July 2017

  Tests the wideband functions with real data derived from input
  speech samples.  Test vectors are dumped to Octave vectors so they can
  be verified against the Octave version of the wideband functions.

  Supports rapid go/no-go testing of the C port when any canges are
  made, and flushes out bugs in the C and octave version.

\*---------------------------------------------------------------------------*/

void test_with_real_data(int argc, char *argv[])
{
    int Nt = C2WB_NT;
    int Fs = C2WB_FS;
    int K = C2WB_K;

    C2CONST c2const = c2const_create(Fs);
    int n_samp = c2const.n_samp;
    int m_pitch = c2const.m_pitch;
    short buf[n_samp]; /* input/output buffer                   */
    float Sn[m_pitch]; /* float input speech samples            */
    COMP Sw[FFT_ENC]; /* DFT of Sn[]                           */
    codec2_fft_cfg fft_fwd_cfg; /* fwd FFT states                        */
    float w[m_pitch]; /* time domain hamming window            */
    COMP W[FFT_ENC]; /* DFT of w[]                            */
    MODEL model[FRAMES];
    void *nlp_states;
    float pitch, prev_f0;
    int i, m, f;

    if (argc != 2) {
        printf("test_with_real_data usage: .%s RawFile\n", argv[0]);
        exit(1);
    }

    nlp_states = nlp_create(&c2const);
    prev_f0 = 1.0 / P_MAX_S;
    fft_fwd_cfg = codec2_fft_alloc(FFT_ENC, 0, NULL, NULL);
    make_analysis_window(&c2const, fft_fwd_cfg, w, W);

    for (i = 0; i < m_pitch; i++) {
        Sn[i] = 1.0;
    }

    WIDEBAND_MAP wb_map;
    setup_map(&wb_map, Nt, K);
    int n_block_frames = C2WB_NT * C2WB_DEC;

    float model_octave[FRAMES][MAX_AMP + 2]; // model params in matrix format, useful for C <-> Octave  
    float rate_K_surface[FRAMES][K]; // rate K vecs for each frame, form a surface that makes pretty graphs
    float rate_K_surface_[FRAMES][K];
    MODEL model_block_[n_block_frames];

    int Nblocks = FRAMES / n_block_frames;
    float dct2_sd[Nblocks];
    int qn;

    for (f = 0; f < FRAMES; f++) {
        for (m = 0; m < MAX_AMP + 2; m++) {
            model_octave[f][m] = 0.0;
        }
    }

    FILE *fin = fopen(argv[1], "rb");
    if (fin == NULL) {
        fprintf(stderr, "Problem opening %s\n", argv[1]);
        exit(1);
    }

    for (f = 0; f < FRAMES; f++) {
        assert(fread(buf, sizeof (short), n_samp, fin) == n_samp);

        /* shift buffer of input samples, and insert new samples */

        for (i = 0; i < m_pitch - n_samp; i++) {
            Sn[i] = Sn[i + n_samp];
        }
        for (i = 0; i < n_samp; i++) {
            Sn[i + m_pitch - n_samp] = buf[i];
        }

        /* Estimate Sinusoidal Model Parameters ----------------------*/

        nlp(nlp_states, Sn, n_samp, &pitch, Sw, W, &prev_f0);
        model[f].Wo = TWO_PI / pitch;

        dft_speech(&c2const, fft_fwd_cfg, Sw, Sn, w);
        two_stage_pitch_refinement(&c2const, &model[f], Sw);
        estimate_amplitudes(&model[f], Sw, W, 1);
        est_voicing_mbe(&c2const, &model[f], Sw, W);

        fprintf(stderr, "f: %d Wo: %4.3f L: %d v: %d\n", f, model[f].Wo, model[f].L, model[f].voiced);

        /* log some vectors for sinusoidal model */

        model_octave[f][0] = model[f].Wo;
        model_octave[f][1] = model[f].L;

        for (m = 1; m <= model[f].L; m++) {
            model_octave[f][m + 1] = model[f].A[m];
        }

        /* once we have collected a block of samples process ----------*/

        if (((f + 1) % n_block_frames) == 0) {

            /* wideband processing ----------------------------------------*/

            int block_st = f - n_block_frames + 1;
            wideband_enc_dec(&c2const, n_block_frames, &model[block_st], &wb_map,
                             model_block_, dct2_sd, &qn,
                             &rate_K_surface[block_st],
                             &rate_K_surface_[block_st]);

            fprintf(stderr, "  Performed a wideband_enc_dec() call, qn: %d\n", qn);

            /* todo: add a synthesis stage here to output decoded speech from model_block_ */
        }
    }

    fclose(fin);

    /* save vectors in Octave format */

    FILE *fout = fopen("tc2wideband_out.txt", "wt");
    assert(fout != NULL);
    fprintf(fout, "# Created by tc2wideband.c\n");
    octave_save_float(fout, "rate_K_surface_c", (float*) rate_K_surface, FRAMES, K, K);
    octave_save_float(fout, "model_c", (float*) model_octave, FRAMES, MAX_AMP + 2, MAX_AMP + 2);
    fclose(fout);

    printf("Done! Now run\n  octave:1> tc2wideband(\"../path/to/tc2wideband_out.txt\")\n");
}
