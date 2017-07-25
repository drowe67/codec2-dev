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

    float y[] = {1, 2, 3, 4};
    float expect_dct[] = {20, -6.30864406, 0, -0.44834153};
    float res_dct[N];

    dct(N, y, res_dct);


    for (i = 0; i < N; i++) {
        //printf("Result = %f : %f\n", res_dct[i].real, res_dct[i].imag);
        //printf("Expected = %f\n", expect_dct[i]);
        if (abs(expect_dct[i] - res_dct[i]) > 0.00001)
            test_failed_f(expect_dct[i], res_dct[i]);
    }

    test("idct");
    float *expect_idct = y;
    float res_idct[N];

    idct(N, res_dct, res_idct);


    for (i = 0; i < N; i++) {
        //printf("Result = %f : %f\n", res_dct[i].real, res_dct[i].imag);
        //printf("Expected = %f\n", expect_dct[i]);
        if (abs(expect_idct[i] - res_idct[i]) > 0.00001)
            test_failed_f(expect_idct[i], res_idct[i]);
    }

    test("dct2");

    int M = 3;
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
    dct2(M, N, y2, res_dct2);



    for (i = 0; i < M; i++) {
        for (j = 0; j < N; j++) {
            //printf("Result = %f : %f\n", res_dct[i].real, res_dct[i].imag);
            //printf("Expected = %f\n", expect_dct[i]);
            if (abs(expect_dct2[i][j] - res_dct2[i][j]) > 0.00001)
                test_failed_f(expect_dct2[i][j], res_dct2[i][j]);
        }
    }
    
    test("idct2");


    float expect_idct2[3][4] = {
        {1, 2, 3, 4},
        {5, 6, 7, 8},
        {3, 1, 2, 3}
    };

    
    float res_idct2[M][N];
    idct2(M, N, res_dct2, res_idct2);

    for (i = 0; i < M; i++) {
        for (j = 0; j < N; j++) {
            //printf("Result = %f : %f\n", res_dct[i].real, res_dct[i].imag);
            //printf("Expected = %f\n", expect_dct[i]);
            if (abs(expect_idct2[i][j] - res_idct2[i][j]) > 0.00001)
                test_failed_f(expect_idct2[i][j], res_idct2[i][j]);
        }
    }
}