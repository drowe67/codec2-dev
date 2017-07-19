#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>

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

float mean(float data[], int n);
int unit_test();
float std(float data[], int rows);
void diff_de(int rows, int cols, float D[rows][cols], float E[rows][cols], float diff_de[rows][cols]);
void array_col_to_row(int rows, int cols, float data[rows][cols], int col, float res[]);
void std_on_cols(int rows, int cols, float data[rows][cols], float res[]);
float mean_std_diff_de(int rows, int cols, float D[rows][cols], float E[rows][cols]);
void test_wideband_enc_dec();
void setup_map(WIDEBAND_MAP * wb_map, int Nt, int K);

char *fn;

void test(char * tfn) {
  fn = tfn;
  printf("test function: %s\n", fn);
}

void test_failed() {
  printf("Failed to calculate %s.\n", fn);
  exit(1);
}

void test_failed_s(char * expected, char * res) {

  printf("Failed to calculate %s.\n", fn);

  printf("expected: %s\ngot: %s\n", expected, res);
  exit(1);
}

void test_failed_f(float expected, float res) {

  printf("Failed to calculate %s.\n", fn);
  printf("expected: %f\ngot: %f\n", expected, res);
  exit(1);
}

int main(int argc, char *argv[]) {

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

  for (int r = 0; r < rows; r++) {
    for (int c = 0; c < cols; c++) {
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
  for (int r = 0; r < rows; r++) {
    for (int c = 0; c < cols; c++) {
      float d = rand();
      data[r][c] = d;
      expect_data[c][r] = d;
    }
  }
  for (int c = 0; c < cols; c++) {

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


  return 1;
}

void test_wideband_enc_dec() {

  test("wideband_enc_dec");

  short *buf;
  unsigned char *bits;
  //int            nsam, nbit, i, r;
  int nbit;
  int nsam;
  int n_block_frames = C2WB_NT * C2WB_DEC;
  int test_me = 123;
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


  for (int n = 0; n < n_block_frames; n++) {
    model_block[n].L = 10;
    model_block[n].Wo = (float) rand() / RAND_MAX;
    for (int i = 0; i < MAX_AMP + 1; i++) {
      model_block[n].phi[i] = (float) rand() / RAND_MAX;
      model_block[n].A[i] = (float) rand() / RAND_MAX;
    }
    model_block[n].voiced = round((float) rand() / RAND_MAX);
  }

  int n = 0;

  printf("setup complete. now calling the function\n");

  wideband_enc_dec(&c2const, n_block_frames, model_block, &wb_map,
          model_block_, dct2_sd, &qn, &rate_K_surface_block[n], &rate_K_surface_block_[n]);

  codec2_destroy(codec2);
  free(bits);
  free(buf);
  printf("made it to the end\n");
}

