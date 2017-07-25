#include "dct2.h"

/*
 * 1-dimensional DCT, taking and returning only the real components of the result.
 * The implementation relies on the real FFT used elsewhere within Codec 2, ensuring
 * that any optimisations of the complex code are reused, rather than reimplemented.
 */
void dct(const int N, float y[N], float res[N])
{

    int i;
    float y2[2 * N]; // input to FFT will be double the size of the input
    COMP c[N + 1];
    COMP phi[N];


    codec2_fftr_cfg cfg = codec2_fftr_alloc(2 * N, 0, NULL, NULL);

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

void dct2(const int M, const int N, float y[M][N], float res[M][N])
{
    float a[M][N];

    int i;
    int j;
    for (i = 0; i < M; i++) {
        dct(N, y[i], a[i]);
    }

    float a_col[M];
    float b_col[M];

    for (j = 0; j < N; j++) {
        for (i = 0; i < M; i++) {
            a_col[i] = a[i][j];
        }

        dct(M, a_col, b_col);

        for (i = 0; i < M; i++) {
            res[i][j] = b_col[i];
        }
    }
}

//def idct(a):
//    N = len(a)
//    c = empty(N+1,complex)
//
//    phi = exp(1j*pi*arange(N)/(2*N))
//    c[:N] = phi*a
//    c[N] = 0.0
//    return irfft(c)[:N]

/* 
 * Inverse DCT 1-dimensional.
 */

void idct(const int N, float a[N], float res[N])
{
    
    int nfft = 2 * N ;
    
    // if nfft, the size of the inverse FFT output is odd, extend it by one
    if(nfft & 1)
        nfft++;
    
    int i;
    COMP ac;
    COMP c[N + 1];
    COMP phi[N];
    float res_fft[nfft];

    codec2_fftr_cfg cfg = codec2_fftr_alloc(nfft, 1, NULL, NULL);

    assert(cfg);
    
    for (i = 0; i < N; i++) {
        float p;
        p = PI * i / (2 * N);
        phi[i] = comp_exp_j(p);
        ac.real = a[i];
        ac.imag = 0;
        c[i] = cmult(phi[i], ac);
    }

    c[N].real = 0;
    c[N].imag = 0;
    
    codec2_fftri(cfg, c, res_fft);
    
    // Scale the result
    for(i=0;i<N; i++){
        res[i] = res_fft[i] / nfft;
    }


}

//def idct2(b):
//    M = b.shape[0]
//    N = b.shape[1]
//    a = empty([M,N],float)
//    y = empty([M,N],float)
//
//    for i in range(M):
//        a[i,:] = idct(b[i,:])
//    for j in range(N):
//        y[:,j] = idct(a[:,j])
//
//    return y

/*
 *  Inverse DCT 2-dimensional
 */
void idct2(int M, int N, float y[M][N], float res[M][N])
{

    float a[M][N];

    int i;
    int j;
    for (i = 0; i < M; i++) {
        idct(N, y[i], a[i]);
    }

    float a_col[M];
    float b_col[M];

    for (j = 0; j < N; j++) {
        for (i = 0; i < M; i++) {
            a_col[i] = a[i][j];
        }

        idct(M, a_col, b_col);

        for (i = 0; i < M; i++) {
            res[i][j] = b_col[i];
        }
    }

}
