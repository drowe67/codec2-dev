/*
  FILE...: mpdecode_core.h
  AUTHOR.: David Rowe
  CREATED: Sep 2016

  C-callable core functions for MpDecode, so they can be used for
  Octave and C programs.  Also some convenience functions to help use
  the C-callable LDPC decoder in C programs.
*/

#ifndef __MPDECODE_CORE__
#define __MPDECODE_CORE__

#include "comp.h"
#include <stdint.h>

struct LDPC {
    int max_iter;
    int dec_type;
    int q_scale_factor;
    int r_scale_factor;
    int CodeLength;
    int NumberParityBits;
    int NumberRowsHcols;
    int max_row_weight;
    int max_col_weight;
    int data_bits_per_frame;
    int coded_bits_per_frame;
    int coded_syms_per_frame;
    uint16_t *H_rows;
    uint16_t *H_cols;
};

extern void ldpc_init(struct LDPC *ldpc, int *size_common);
extern void ldpc_free_mem(struct LDPC *ldpc);

extern void encode(struct LDPC *ldpc, unsigned char ibits[], unsigned char pbits[]);

int run_ldpc_decoder(struct LDPC *ldpc, void *ldoc_common, char out_char[], 
			float input[], int *parityCheckCount);
extern void ldpc_dump_nodes(struct LDPC *ldpc);

extern void sd_to_llr(float llr[], double sd[], int n);

extern void Demod2D(double symbol_likelihood[], COMP r[], COMP S_matrix[], float EsNo, float fading[], float mean_amp, int number_symbols);
extern void Somap(double bit_likelihood[], double symbol_likelihood[], int number_symbols);
extern void symbols_to_llrs(double llr[], COMP rx_qpsk_symbols[], float rx_amps[], float EsNo, float mean_amp, int nsyms);

#endif
