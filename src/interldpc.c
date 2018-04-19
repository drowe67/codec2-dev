/*---------------------------------------------------------------------------*\

  FILE........: interldpc.c
  AUTHOR......: David Rowe
  DATE CREATED: April 2018

  Helper functions for interleaved LDPC modems.

\*---------------------------------------------------------------------------*/

/*
  Copyright (C) 2018 David Rowe

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

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "interldpc.h"
#include "codec2_ofdm.h"
#include "ofdm_internal.h"
#include "mpdecode_core.h"
#include "gp_interleaver.h"
#include "HRA_112_112.h"

/* CRC type function, used to compare QPSK vectors when debugging */

COMP test_acc(COMP v[], int n) {
    COMP acc = {0.0,0.0};
    int i;
    for(i=0; i<n; i++) {
        acc.real += round(v[i].real);
        acc.imag += round(v[i].imag);
        //fprintf(stderr, "%d %10f %10f %10f %10f\n", i, round(v[i].real), round(v[i].imag), acc.real, acc.imag);
    }
    return acc;
}

void printf_n(COMP v[], int n) {
    int i;
    for(i=0; i<n; i++) {
        fprintf(stderr, "%d %10f %10f\n", i, round(v[i].real), round(v[i].imag));
    }
}

void set_up_hra_112_112(struct LDPC *ldpc) {
    ldpc->max_iter = HRA_112_112_MAX_ITER;
    ldpc->dec_type = 0;
    ldpc->q_scale_factor = 1;
    ldpc->r_scale_factor = 1;
    ldpc->CodeLength = HRA_112_112_CODELENGTH;
    ldpc->NumberParityBits = HRA_112_112_NUMBERPARITYBITS;
    ldpc->NumberRowsHcols = HRA_112_112_NUMBERROWSHCOLS;
    ldpc->max_row_weight = HRA_112_112_MAX_ROW_WEIGHT;
    ldpc->max_col_weight = HRA_112_112_MAX_COL_WEIGHT;
    ldpc->H_rows = HRA_112_112_H_rows;
    ldpc->H_cols = HRA_112_112_H_cols;

    /* provided for convenience and to match Octave vaiable names */
    
    ldpc->data_bits_per_frame = HRA_112_112_CODELENGTH - HRA_112_112_NUMBERPARITYBITS;
    ldpc->coded_bits_per_frame = HRA_112_112_CODELENGTH;
    ldpc->coded_syms_per_frame = ldpc->coded_bits_per_frame/OFDM_BPS;
}

void ldpc_encode_frame(struct LDPC *ldpc, int codeword[], unsigned char tx_bits_char[]) {
    unsigned char pbits[ldpc->NumberParityBits];
    int           i,j;
    
    encode(ldpc, tx_bits_char, pbits);
    for(i=0; i<ldpc->data_bits_per_frame; i++) {
        codeword[i] = tx_bits_char[i];
    }
    for(j=0; i<ldpc->coded_bits_per_frame; i++,j++) {
        codeword[i] = pbits[j];
    }
}

void qpsk_modulate_frame(COMP tx_symbols[], int codeword[], int n) {
    int s,i;
    int dibit[2];
    complex float qpsk_symb;
    
    for(s=0,i=0; i<n; s += 2,i++) {
        dibit[0] = codeword[s+1] & 0x1;
        dibit[1] = codeword[s] & 0x1;
        qpsk_symb = qpsk_mod(dibit);
        tx_symbols[i].real = crealf(qpsk_symb);
        tx_symbols[i].imag = cimagf(qpsk_symb);
    }    
}
