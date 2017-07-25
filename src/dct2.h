/* 
 * File:   dct2.h
 * Author: Phil Ayres
 *
 * Created on 21 July 2017, 14:24
 */

#ifndef DCT2_H
#define	DCT2_H

#include "codec2_fft.h"
#include "comp.h"
#include "comp_prim.h"

void dct(const int N, float y[], float res[]);
void dct2(const int M, const int N, float y[M][N], float res[M][N]);
void idct(const int N, float a[N], float res[N]);
void idct2(int M, int N, float y[M][N], float res[M][N]);

#endif	/* DCT2_H */

