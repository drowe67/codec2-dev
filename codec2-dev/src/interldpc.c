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
#include <string.h>
#include <math.h>

#include "interldpc.h"
#include "codec2_ofdm.h"
#include "mpdecode_core.h"
#include "gp_interleaver.h"
#include "HRA_112_112.h"
#include "test_bits_ofdm.h"

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

void interleaver_sync_state_machine(struct OFDM *ofdm,
                                    struct LDPC *ldpc,
                                    COMP codeword_symbols_de[],
                                    float codeword_amps_de[],
                                    float EsNo, int interleave_frames,
                                    int *iter, int *parityCheckCount, int *Nerrs_coded)
{
    int coded_syms_per_frame = ldpc->coded_syms_per_frame;
    int coded_bits_per_frame = ldpc->coded_bits_per_frame;
    int data_bits_per_frame = ldpc->data_bits_per_frame;
    double llr[coded_bits_per_frame];
    char out_char[coded_bits_per_frame];                    
    char next_sync_state_interleaver[OFDM_STATE_STR];
    
    strcpy(next_sync_state_interleaver, ofdm->sync_state_interleaver);
    if ((strcmp(ofdm->sync_state_interleaver,"search") == 0) && (ofdm->frame_count >= (interleave_frames-1))) {
        symbols_to_llrs(llr, codeword_symbols_de, codeword_amps_de, EsNo, ofdm->mean_amp, coded_syms_per_frame);               
        iter[0] =  run_ldpc_decoder(ldpc, out_char, llr, parityCheckCount);
        Nerrs_coded[0] = data_bits_per_frame - parityCheckCount[0];
        //for(i=0; i<20; i++)
        //    fprintf(stderr,"%d ", out_char[i]);
        //fprintf(stderr,"\n");
        //fprintf(stderr, "     iter: %d pcc: %d Nerrs: %d\n", iter[0], parityCheckCount[0], Nerrs_coded[0]);
        if ((Nerrs_coded[0] == 0) || (interleave_frames == 1)) {
            /* sucessful decode! */
            strcpy(next_sync_state_interleaver, "synced");
            ofdm->frame_count_interleaver = interleave_frames;
        }
    }
    strcpy(ofdm->sync_state_interleaver, next_sync_state_interleaver);
}


/* measure uncoded (raw) bit errors over interleaver frame */

int count_uncoded_errors(struct LDPC *ldpc, int Nerrs_raw[], int interleave_frames, COMP codeword_symbols_de[])
{
    int i,j,Nerrs,Terrs;

    int coded_syms_per_frame = ldpc->coded_syms_per_frame;
    int coded_bits_per_frame = ldpc->coded_bits_per_frame;
    int rx_bits_raw[coded_bits_per_frame];

    assert(sizeof(test_codeword)/sizeof(int) == coded_bits_per_frame);

    Terrs = 0;
    for (j=0; j<interleave_frames; j++) {
        for(i=0; i<coded_syms_per_frame; i++) {
            int bits[2];
            complex float s = codeword_symbols_de[j*coded_syms_per_frame+i].real + I*codeword_symbols_de[j*coded_syms_per_frame+i].imag;
            qpsk_demod(s, bits);
            rx_bits_raw[OFDM_BPS*i]   = bits[1];
            rx_bits_raw[OFDM_BPS*i+1] = bits[0];
        }
        Nerrs = 0;
        for(i=0; i<coded_bits_per_frame; i++) {
            //fprintf(stderr, "%d %d %d\n", i, test_codeword[i], rx_bits_raw[i]);
            if (test_codeword[i] != rx_bits_raw[i]) {
                Nerrs++;
            }
        }
                                
        Nerrs_raw[j] = Nerrs;
        Terrs += Nerrs;
    }

    return Terrs;
}

int count_errors(int tx_bits[], char rx_bits[], int n) {
    int i;
    int Nerrs = 0;
    
    Nerrs = 0;
    for(i=0; i<n; i++) {
        if (tx_bits[i] != rx_bits[i]) {
            Nerrs++;
        }
    }
    return Nerrs;
}


#ifdef NOT_USED
/* UW never changes so we can pre-load tx_symbols with modulated UW */

void build_modulated_uw(struct OFDM *ofdm, complex float tx_symbols[], uint8_t txt_bits[])
{
    int  uw_txt_bits[OFDM_NUWBITS+OFDM_NTXTBITS];
    COMP uw_txt_syms[(OFDM_NUWBITS+OFDM_NTXTBITS)/OFDM_BPS];
    int i,j;

    for(i=0; i<OFDM_NUWBITS; i++) {
        uw_txt_bits[i] = ofdm->tx_uw[i];
    }
    /* clear txt bits for now, they can be added in later */
    for(j=0; j<OFDM_NTXTBITS; j++,i++) {
        uw_txt_bits[i] = txt_bits[j];
        //fprintf(stderr, "txt_bits[%d] = %d\n", j, txt_bits[j]);
    }    
    qpsk_modulate_frame(uw_txt_syms, uw_txt_bits, (OFDM_NUWBITS+OFDM_NTXTBITS)/OFDM_BPS);
    for(i=0; i<(OFDM_NUWBITS+OFDM_NTXTBITS)/OFDM_BPS; i++) {
        tx_symbols[i] = uw_txt_syms[i].real + I * uw_txt_syms[i].imag;
    }
}
#endif

/* 
   Given an array of tx_bits, LDPC encodes, interleaves, and OFDM
   modulates.

   Note this could be refactored to save memory, e.g. for embedded
   applications we could call ofdm_txframe on a frame by frame
   basis 
*/

void ofdm_ldpc_interleave_tx(struct OFDM *ofdm, struct LDPC *ldpc, complex float tx_sams[], uint8_t tx_bits[], uint8_t txt_bits[], int interleave_frames)
{
    int coded_syms_per_frame = ldpc->coded_syms_per_frame;
    int coded_bits_per_frame = ldpc->coded_bits_per_frame;
    int data_bits_per_frame = ldpc->data_bits_per_frame;

    int codeword[coded_bits_per_frame];
    COMP coded_symbols[interleave_frames*coded_syms_per_frame];
    COMP coded_symbols_inter[interleave_frames*coded_syms_per_frame];
    int Nsamperframe = ofdm_get_samples_per_frame();
    complex float tx_symbols[OFDM_BITSPERFRAME/OFDM_BPS];
    int j;
    
    for (j=0; j<interleave_frames; j++) {
        ldpc_encode_frame(ldpc, codeword, &tx_bits[j*data_bits_per_frame]);
        qpsk_modulate_frame(&coded_symbols[j*coded_syms_per_frame], codeword, coded_syms_per_frame);
    }
    gp_interleave_comp(coded_symbols_inter, coded_symbols, interleave_frames*coded_syms_per_frame);
    for (j=0; j<interleave_frames; j++) {            
        ofdm_assemble_modem_frame(tx_symbols, &coded_symbols_inter[j*coded_syms_per_frame], &txt_bits[OFDM_NTXTBITS*j]);
        ofdm_txframe(ofdm, &tx_sams[j*Nsamperframe], tx_symbols);
    }
}


