/*---------------------------------------------------------------------------*\

  FILE........: interldpc.h
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

#ifndef __INTERLDPC__
#define __INTERLDPC__

#include "comp.h"
#include "mpdecode_core.h"
#include "ofdm_internal.h"

/* CRC type function, used to compare QPSK vectors when debugging */

COMP test_acc(COMP v[], int n);
void printf_n(COMP v[], int n);
void set_up_hra_112_112(struct LDPC *ldpc);
void ldpc_encode_frame(struct LDPC *ldpc, int codeword[], unsigned char tx_bits_char[]);
void qpsk_modulate_frame(COMP tx_symbols[], int codeword[], int n);
void interleaver_sync_state_machine(struct OFDM *ofdm,
                                    struct LDPC *ldpc,
                                    COMP codeword_symbols_de[],
                                    float codeword_amps_de[],
                                    float EsNo, int interleave_frames,
                                    int *inter, int *parityCheckCount, int *Nerrs_coded);
int count_uncoded_errors(struct LDPC *ldpc, int Nerrs_raw[], int interleave_frames, COMP codeword_symbols_de[]);

#endif
