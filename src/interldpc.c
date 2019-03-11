/*---------------------------------------------------------------------------*\

  FILE........: interldpc.c
  AUTHOR......: David Rowe
  DATE CREATED: April 2018

  Helper functions for interleaved LDPC waveforms.

\*---------------------------------------------------------------------------*/

/*
  Copyright (C) 2018 David Rowe

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
#include <stdint.h>
#include <string.h>
#include <math.h>

#include "interldpc.h"
#include "ofdm_internal.h"
#include "mpdecode_core.h"
#include "gp_interleaver.h"
#include "HRA_112_112.h"
#include "HRAb_396_504.h"

/* CRC type function, used to compare QPSK vectors when debugging */

COMP test_acc(COMP v[], int n) {
    COMP acc = {0.0f, 0.0f};
    int i;
    
    for (i = 0; i < n; i++) {
        acc.real += roundf(v[i].real);
        acc.imag += roundf(v[i].imag);
    }

    return acc;
}

void printf_n(COMP v[], int n) {
    int i;
    
    for (i = 0; i < n; i++) {
        fprintf(stderr, "%d %10f %10f\n", i, round(v[i].real), round(v[i].imag));
    }
}

void set_up_ldpc_constants(struct LDPC *ldpc, int code_length, int parity_bits, int bps) {
    /* following provided for convenience and to match Octave variable names */

    /* these remain fixed */
    ldpc->ldpc_data_bits_per_frame = code_length - parity_bits;
    ldpc->ldpc_coded_bits_per_frame = code_length;

    /* in the case there are some unused data bits, these may be
       modified to be less that ldpc->ldpc_xxx versions above. We
       place known bits in the unused data bit positions, which make
       the code stronger, and allow us to mess with different speech
       codec bit allocations without designing new LDPC codes. */
    
    ldpc->data_bits_per_frame = ldpc->ldpc_data_bits_per_frame;
    ldpc->coded_bits_per_frame = ldpc->ldpc_coded_bits_per_frame;
    ldpc->coded_syms_per_frame = ldpc->coded_bits_per_frame / bps;
}

void set_data_bits_per_frame(struct LDPC *ldpc, int new_data_bits_per_frame, int bps) {
    ldpc->data_bits_per_frame = new_data_bits_per_frame;
    ldpc->coded_bits_per_frame = ldpc->data_bits_per_frame + ldpc->NumberParityBits;
    ldpc->coded_syms_per_frame = ldpc->coded_bits_per_frame / bps;
}
    
// TODO: this should be in (n,k) = (224,112) format, fix some time

void set_up_hra_112_112(struct LDPC *ldpc, struct OFDM_CONFIG *config) {
    ldpc->max_iter = HRA_112_112_MAX_ITER;
    ldpc->dec_type = 0;
    ldpc->q_scale_factor = 1;
    ldpc->r_scale_factor = 1;
    ldpc->CodeLength = HRA_112_112_CODELENGTH;
    ldpc->NumberParityBits = HRA_112_112_NUMBERPARITYBITS;
    ldpc->NumberRowsHcols = HRA_112_112_NUMBERROWSHCOLS;
    ldpc->max_row_weight = HRA_112_112_MAX_ROW_WEIGHT;
    ldpc->max_col_weight = HRA_112_112_MAX_COL_WEIGHT;
    ldpc->H_rows = (uint16_t *) HRA_112_112_H_rows;
    ldpc->H_cols = (uint16_t *) HRA_112_112_H_cols;

    /* provided for convenience and to match Octave variable names */

    set_up_ldpc_constants(ldpc, HRA_112_112_CODELENGTH, HRA_112_112_NUMBERPARITYBITS, config->bps);
}

// Note code #defines below should be in (n,k) = (504,396)
// TODO : fix this some time
void set_up_hra_504_396(struct LDPC *ldpc, struct OFDM_CONFIG *config) {
    ldpc->max_iter = HRAb_396_504_MAX_ITER;
    ldpc->dec_type = 0;
    ldpc->q_scale_factor = 1;
    ldpc->r_scale_factor = 1;
    ldpc->CodeLength = HRAb_396_504_CODELENGTH;
    ldpc->NumberParityBits = HRAb_396_504_NUMBERPARITYBITS;
    ldpc->NumberRowsHcols = HRAb_396_504_NUMBERROWSHCOLS;
    ldpc->max_row_weight = HRAb_396_504_MAX_ROW_WEIGHT;
    ldpc->max_col_weight = HRAb_396_504_MAX_COL_WEIGHT;
    ldpc->H_rows = (uint16_t *) HRAb_396_504_H_rows;
    ldpc->H_cols = (uint16_t *) HRAb_396_504_H_cols;

    set_up_ldpc_constants(ldpc, HRAb_396_504_CODELENGTH, HRAb_396_504_NUMBERPARITYBITS, config->bps);
}

void ldpc_encode_frame(struct LDPC *ldpc, int codeword[], unsigned char tx_bits_char[]) {
    unsigned char pbits[ldpc->NumberParityBits];
    int i, j;

    if (ldpc->data_bits_per_frame == ldpc->ldpc_data_bits_per_frame) {
        /* we have enough data bits to fill the codeword */
        encode(ldpc, tx_bits_char, pbits);
    } else {        
        unsigned char tx_bits_char_padded[ldpc->ldpc_data_bits_per_frame];
        /* some unused data bits, set these to known values to strengthen code */    
        memcpy(tx_bits_char_padded, tx_bits_char, ldpc->data_bits_per_frame);
        for (i = ldpc->data_bits_per_frame; i < ldpc->ldpc_data_bits_per_frame; i++)
            tx_bits_char_padded[i] = 1;
        encode(ldpc, tx_bits_char_padded, pbits);
    }
          
    /* output codeword is concatenation of (used) data bits and parity
       bits, we don't bother sending unused (known) data bits */
    for (i = 0; i < ldpc->data_bits_per_frame; i++) {
        codeword[i] = tx_bits_char[i];
    }
    for (j = 0; j < ldpc->NumberParityBits; i++, j++) {
        codeword[i] = pbits[j];
    }
}

void qpsk_modulate_frame(COMP tx_symbols[], int codeword[], int n) {
    int s, i;
    int dibit[2];
    complex float qpsk_symb;

    for (s = 0, i = 0; i < n; s += 2, i++) {
        dibit[0] = codeword[s + 1] & 0x1;
        dibit[1] = codeword[s] & 0x1;
        qpsk_symb = qpsk_mod(dibit);
        tx_symbols[i].real = crealf(qpsk_symb);
        tx_symbols[i].imag = cimagf(qpsk_symb);
    }
}

void interleaver_sync_state_machine(struct OFDM *ofdm,
        struct LDPC *ldpc,
        struct OFDM_CONFIG *config,
        COMP codeword_symbols_de[],
        float codeword_amps_de[],
        float EsNo, int interleave_frames,
        int *iter, int *parityCheckCount, int *Nerrs_coded) {
    int coded_syms_per_frame = ldpc->coded_syms_per_frame;
    int coded_bits_per_frame = ldpc->coded_bits_per_frame;
    int data_bits_per_frame = ldpc->data_bits_per_frame;
    float llr[coded_bits_per_frame];
    uint8_t out_char[coded_bits_per_frame];
    State next_sync_state_interleaver;

    next_sync_state_interleaver = ofdm->sync_state_interleaver;

    if ((ofdm->sync_state_interleaver == search) && (ofdm->frame_count >= (interleave_frames - 1))) {
        symbols_to_llrs(llr, codeword_symbols_de, codeword_amps_de, EsNo, ofdm->mean_amp, coded_syms_per_frame);
        iter[0] = run_ldpc_decoder(ldpc, out_char, llr, parityCheckCount);
        Nerrs_coded[0] = data_bits_per_frame - parityCheckCount[0];

        if ((Nerrs_coded[0] == 0) || (interleave_frames == 1)) {
            /* sucessful decode! */
            next_sync_state_interleaver = synced;
            ofdm->frame_count_interleaver = interleave_frames;
        }
    }

    ofdm->sync_state_interleaver = next_sync_state_interleaver;
}

/* measure uncoded (raw) bit errors over interleaver frame, note we
   don't include txt bits as this is done after we dissassemmble the
   frame */

int count_uncoded_errors(struct LDPC *ldpc, struct OFDM_CONFIG *config, int Nerrs_raw[], int interleave_frames, COMP codeword_symbols_de[]) {
    int i, j, Nerrs, Terrs;

    int coded_syms_per_frame = ldpc->coded_syms_per_frame;
    int coded_bits_per_frame = ldpc->coded_bits_per_frame;
    int data_bits_per_frame = ldpc->data_bits_per_frame;
    int rx_bits_raw[coded_bits_per_frame];

    /* generate test codeword from known payload data bits */

    int test_codeword[coded_bits_per_frame];
    uint16_t r[data_bits_per_frame];
    uint8_t tx_bits[data_bits_per_frame];

    ofdm_rand(r, data_bits_per_frame);
    
    for (i = 0; i < data_bits_per_frame; i++) {
        tx_bits[i] = r[i] > 16384;
    }
    
    ldpc_encode_frame(ldpc, test_codeword, tx_bits);

    Terrs = 0;
    for (j = 0; j < interleave_frames; j++) {
        for (i = 0; i < coded_syms_per_frame; i++) {
            int bits[2];
            complex float s = codeword_symbols_de[j * coded_syms_per_frame + i].real + I * codeword_symbols_de[j * coded_syms_per_frame + i].imag;
            qpsk_demod(s, bits);
            rx_bits_raw[config->bps * i] = bits[1];
            rx_bits_raw[config->bps * i + 1] = bits[0];
        }
        
        Nerrs = 0;
        
        for (i = 0; i < coded_bits_per_frame; i++) {
            if (test_codeword[i] != rx_bits_raw[i]) {
                Nerrs++;
            }
        }

        Nerrs_raw[j] = Nerrs;
        Terrs += Nerrs;
    }

    return Terrs;
}

int count_errors(uint8_t tx_bits[], uint8_t rx_bits[], int n) {
    int i;
    int Nerrs = 0;

    for (i = 0; i < n; i++) {
        if (tx_bits[i] != rx_bits[i]) {
            Nerrs++;
        }
    }
    
    return Nerrs;
}

/*
   Given an array of tx_bits, LDPC encodes, interleaves, and OFDM
   modulates.

   Note this could be refactored to save memory, e.g. for embedded
   applications we could call ofdm_txframe on a frame by frame
   basis
 */

void ofdm_ldpc_interleave_tx(struct OFDM *ofdm, struct LDPC *ldpc, complex float tx_sams[], uint8_t tx_bits[], uint8_t txt_bits[], int interleave_frames, struct OFDM_CONFIG *config) {
    int coded_syms_per_frame = ldpc->coded_syms_per_frame;
    int coded_bits_per_frame = ldpc->coded_bits_per_frame;
    int data_bits_per_frame = ldpc->data_bits_per_frame;
    int ofdm_bitsperframe = ofdm_get_bits_per_frame();
    int codeword[coded_bits_per_frame];
    COMP coded_symbols[interleave_frames * coded_syms_per_frame];
    COMP coded_symbols_inter[interleave_frames * coded_syms_per_frame];
    int Nsamperframe = ofdm_get_samples_per_frame();
    complex float tx_symbols[ofdm_bitsperframe / config->bps];
    int j;

    for (j = 0; j < interleave_frames; j++) {
        ldpc_encode_frame(ldpc, codeword, &tx_bits[j * data_bits_per_frame]);
        qpsk_modulate_frame(&coded_symbols[j * coded_syms_per_frame], codeword, coded_syms_per_frame);
    }

    gp_interleave_comp(coded_symbols_inter, coded_symbols, interleave_frames * coded_syms_per_frame);

    for (j = 0; j < interleave_frames; j++) {
        ofdm_assemble_modem_frame_symbols(tx_symbols, &coded_symbols_inter[j * coded_syms_per_frame], &txt_bits[config->txtbits * j]);
        ofdm_txframe(ofdm, &tx_sams[j * Nsamperframe], tx_symbols);
    }
}

