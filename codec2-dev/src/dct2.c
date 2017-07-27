/*---------------------------------------------------------------------------*\

  FILE........: dct2.c
  AUTHOR......: Phil Ayres
  DATE CREATED: July 2017

 * DCT functions based on existing Codec 2 FFT
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

#include "dct2.h"


/*
 * Call dct_config with the number of parameters to be supplied to a subsequent 
 * DCT or DCT2 function. Returns a configuration that must be
 * passed as the initial argument to the subsequent calls.
 */
codec2_dct_cfg dct_config(int P)
{
    return codec2_fftr_alloc(P * 2, 0, NULL, NULL);
}

/*
 * Call idct_config with the number of parameters to be supplied to a subsequent 
 * inverse DCT or inverse DCT2 function. Returns a configuration that must be
 * passed as the initial argument to the subsequent calls.
 */
codec2_dct_cfg idct_config(int P)
{
    return codec2_fftr_alloc(P * 2, 1, NULL, NULL);
}

/*
 * Free memory from previous dct_config() and idct_config() calls
 */
void dct_cfg_free(codec2_dct_cfg cfg)
{
   codec2_fft_free((codec2_fft_cfg)cfg);
}


/*
 * 1-dimensional DCT, taking and returning only the real components of the result.
 * The implementation relies on the real FFT used elsewhere within Codec 2, ensuring
 * that any optimisations of the complex code are reused, rather than reimplemented.
 */
void dct(codec2_dct_cfg cfg, const int N, float y[N], float res[N])
{

    int i;
    float y2[2 * N]; // input to FFT will be double the size of the input
    COMP c[N + 1];
    COMP phi[N];

    for (i = 0; i < N; i++) {
        y2[i] = y[i];
        y2[N + i] = y[N - 1 - i];
    }

    codec2_fftr(cfg, y2, c);

    for (i = 0; i < N; i++) {
        float p;
        if (i == 0)
            p = 0;
        else
            p = -1 * PI * i / (2 * N);
        phi[i] = comp_exp_j(p);
        res[i] = cmult(phi[i], c[i]).real;
    }

}

/* 
 * 2-dimensional DCT, using and returning only the real components.
 * Built on the 1-D DCT.
 */

void dct2(codec2_dct_cfg cfg_m, codec2_dct_cfg cfg_n, const int M, const int N, float y[M][N], float res[M][N])
{
    float a[M][N];

    int i;
    int j;
    for (i = 0; i < M; i++) {
        dct(cfg_n, N, y[i], a[i]);
    }

    
    float a_tran[N][M];

    float b_tran[N][M];
//    for (j = 0; j < N; j++) {
//        for (i = 0; i < M; i++) {
//            a_col[i] = a[i][j];
//        }
//
//        dct(cfg, M, a_col, b_col);
//
//        for (i = 0; i < M; i++) {
//            res[i][j] = b_col[i];
//        }
//    }
    
//    for i in range(M):
  //      a[i,:] = dct(y[i,:])
    
    for(i=0; i<M; i++){        
        dct(cfg_n, N, y[i], a[i]);        
        
        for(j=0; j<N;j++){
            a_tran[j][i] = a[i][j];
        }
        
//        printf("y[i,:]= ");
//        for(j=0; j<N;j++){
//            printf("%f ", y[i][j]);
//        }
//        printf("\n");
//        
//        printf("a[i,:]= ");
//        for(j=0; j<N;j++){
//            printf("%f ", a[i][j]);
//        }
//        printf("\n");
        
    }
    
    
    
//    for j in range(N):
  //      b[:,j] = dct(a[:,j])
    

    for(j=0; j<N;j++){
        
        dct(cfg_m, M, a_tran[j], b_tran[j]);
        
        for(i=0; i<M; i++){
            res[i][j] = b_tran[j][i];
        }
        
//        printf("a[:,j]= ");
//        for(i=0; i<M;i++){
//            printf("%f ", a_tran[j][i]);
//        }
//        printf("\n");
//        
//        printf("b[:,j]= ");
//        for(i=0; i<M;i++){
//            printf("%f ", b_tran[j][i]);
//        }
//        printf("\n");
        
    }
    
    
}


/* 
 * Inverse DCT 1-dimensional.
 */

void idct(codec2_dct_cfg cfg, const int N, float a[N], float res[N])
{
    
    int nfft = N * 2;
    int i;
    COMP ac;
    COMP c[N + 1];
    COMP phi[N];
    float res_fft[nfft];

    assert(cfg);
    
    for (i = 0; i < N; i++) {
        float p;
        p = PI * i / (2 * N);
        phi[i] = comp_exp_j(p);
        //printf("phi: %f+%fj\n", phi[i].real, phi[i].imag);
        ac.real = a[i];
        ac.imag = 0;
        c[i] = cmult(phi[i], ac);
        //printf("c: %f+%fj\n", c[i].real, c[i].imag);
    }

    c[N].real = 0;
    c[N].imag = 0;
    
    codec2_fftri(cfg, c, res_fft);
    
    // Scale the result
    for(i=0;i<N; i++){
        res[i] = res_fft[i] / nfft;
        //printf("res: %f\n", res[i]);
    }


}

/*
 *  Inverse DCT 2-dimensional
 */
void idct2(codec2_dct_cfg cfg_m, codec2_dct_cfg cfg_n, int M, int N, float y[M][N], float res[M][N])
{

    float a[M][N];

    int i;
    int j;
    for (i = 0; i < M; i++) {
        idct(cfg_n, N, y[i], a[i]);
    }

    float a_col[M];
    float b_col[M];

    for (j = 0; j < N; j++) {
        for (i = 0; i < M; i++) {
            a_col[i] = a[i][j];
        }

        idct(cfg_m, M, a_col, b_col);

        for (i = 0; i < M; i++) {
            res[i][j] = b_col[i];
        }
    }

}
