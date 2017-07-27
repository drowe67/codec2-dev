/*---------------------------------------------------------------------------*\

  FILE........: tdct2.c
  AUTHOR......: Phil Ayres
  DATE CREATED: July 2017

 * Unit test for DCT & DCT2 functions
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

#include "t_helpers.h"
#include "dct2.h"

int main()
{
    int i;
    int j;

    test("dct");

    int N = 4;

    codec2_dct_cfg dct_cfg = dct_config(N);

    float y[] = {1, 2, 3, 4};
    float expect_dct[] = {20, -6.30864406, 0, -0.44834153};
    float res_dct[N];

    dct(dct_cfg, N, y, res_dct);


    for (i = 0; i < N; i++) {
        //printf("Result = %f : %f\n", res_dct[i].real, res_dct[i].imag);
        //printf("Expected = %f\n", expect_dct[i]);
        if (abs(expect_dct[i] - res_dct[i]) > 0.00001)
            test_failed_f(expect_dct[i], res_dct[i]);
    }

    dct_cfg_free(dct_cfg);

    test("idct");

    codec2_dct_cfg idct_cfg = idct_config(N);

    float *expect_idct = y;
    float res_idct[N];

    idct(idct_cfg, N, res_dct, res_idct);


    for (i = 0; i < N; i++) {
        //printf("Result = %f : %f\n", res_dct[i].real, res_dct[i].imag);
        //printf("Expected = %f\n", expect_dct[i]);
        if (abs(expect_idct[i] - res_idct[i]) > 0.00001)
            test_failed_f(expect_idct[i], res_idct[i]);
    }


    dct_cfg_free(idct_cfg);

    test("dct2");

    int M = 3;

    codec2_dct_cfg dct_cfg_n = dct_config(N);
    codec2_dct_cfg dct_cfg_m = dct_config(M);


    float expect_dct2[3][4] = {
        { 180, -26.76530997, 8.48528137, 1.90215201},
        { 3.46410162, -9.60123774, -7.34846923, -3.97696289},
        { -66, 5.5432772, 4.24264069, 2.29610059}
    };

    float y2[3][4] = {
        {1, 2, 3, 4},
        {5, 6, 7, 8},
        {3, 1, 2, 3}
    };


    float res_dct2[M][N];
    dct2(dct_cfg_m, dct_cfg_n, M, N, y2, res_dct2);

    /*
        printf("Result\n");
        for (i = 0; i < M; i++) {
            for (j = 0; j < N; j++) {
                printf("%f ", res_dct2[i][j]);
            
                //if (abs(expect_dct2[i][j] - res_dct2[i][j]) > 0.00001)
                  //  test_failed_f(expect_dct2[i][j], res_dct2[i][j]);
            }
            printf("\n");
        }
        printf("Expected\n");
        for (i = 0; i < M; i++) {
            for (j = 0; j < N; j++) {
                printf("%f ", expect_dct2[i][j]);
            
                //if (abs(expect_dct2[i][j] - res_dct2[i][j]) > 0.00001)
                  //  test_failed_f(expect_dct2[i][j], res_dct2[i][j]);
            }
            printf("\n");
        }
     */
    for (i = 0; i < M; i++) {
        for (j = 0; j < N; j++) {

            if (abs(expect_dct2[i][j] - res_dct2[i][j]) > 0.00001)
                test_failed_f(expect_dct2[i][j], res_dct2[i][j]);
        }

    }



    dct_cfg_free(dct_cfg_m);
    dct_cfg_free(dct_cfg_n);
    test("idct2");

    codec2_dct_cfg idct_cfg_n = idct_config(N);
    codec2_dct_cfg idct_cfg_m = idct_config(M);


    float expect_idct2[3][4] = {
        {1, 2, 3, 4},
        {5, 6, 7, 8},
        {3, 1, 2, 3}
    };


    float res_idct2[M][N];
    idct2(idct_cfg_m, idct_cfg_n, M, N, res_dct2, res_idct2);

    for (i = 0; i < M; i++) {
        for (j = 0; j < N; j++) {
            //printf("Result = %f \n", res_idct2[i][j]);
            //printf("Expected = %f\n", expect_idct2[i][j]);
            if (abs(expect_idct2[i][j] - res_idct2[i][j]) > 0.00001)
                test_failed_f(expect_idct2[i][j], res_idct2[i][j]);
        }
    }

    dct_cfg_free(idct_cfg_m);
    dct_cfg_free(idct_cfg_n);


    printf("OK!\n");
}