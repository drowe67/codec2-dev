/*---------------------------------------------------------------------------*\
                                                                             
  FILE........: quantise.h
  AUTHOR......: David Rowe                                                          
  DATE CREATED: 31/5/92                                                       
                                                                             
  Quantisation functions for the sinusoidal coder.  
                                                                             
\*---------------------------------------------------------------------------*/

/*
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

#ifndef __QUANTISE__
#define __QUANTISE__

#define WO_BITS     7
#define WO_LEVELS   (1<<WO_BITS)
#define WO_DT_BITS  3

#define E_BITS      5
#define E_LEVELS    (1<<E_BITS)
#define E_MIN_DB   -10.0
#define E_MAX_DB    40.0

#define LSP_SCALAR_INDEXES    10
#define LSP_DIFF_FREQ_INDEXES  5
#define LSP_DIFF_TIME_BITS     7

#define LSPDT_ALL   0
#define LSPDT_LOW   1
#define LSPDT_HIGH  2

void quantise_init();
float lpc_model_amplitudes(float Sn[], float w[], MODEL *model, int order,
			   int lsp,float ak[]);
void aks_to_M2(float ak[], int order, MODEL *model, float E, float *snr, 
	       int dump);

int   encode_Wo(float Wo);
float decode_Wo(int index);
int   encode_Wo_dt(float Wo, float prev_Wo);
float decode_Wo_dt(int index, float prev_Wo);
void  encode_lsps_scalar(int indexes[], float lsp[], int order);
void  decode_lsps_scalar(float lsp[], int indexes[], int order);
void  encode_lsps_diff_freq_vq(int indexes[], float lsp[], int order);
void  decode_lsps_diff_freq_vq(float lsp_[], int indexes[], int order);
void encode_lsps_diff_time_vq(int indexes[], 
			      float lsp[], 
			      float lsp__prev[], 
			      int order);
void decode_lsps_diff_time_vq(
			      float lsp_[], 
			      int indexes[], 
			      float lsp__prev[],
			      int order);

void lspd_quantise(float lsp[], float lsp_[], int order);
void lspvq_quantise(float lsp[], float lsp_[], int order); 
void lspjnd_quantise(float lsp[], float lsp_[], int order);
void lspdt_quantise(float lsps[], float lsps_[], float lsps__prev[], int mode);

int encode_energy(float e);
float decode_energy(int index);

void pack(unsigned char * bits, unsigned int *nbit, int index, unsigned int index_bits);
int  unpack(const unsigned char * bits, unsigned int *nbit, unsigned int index_bits);

int lsp_bits(int i);
int lspd_bits(int i);

void apply_lpc_correction(MODEL *model);
float speech_to_uq_lsps(float lsp[],
			float ak[],
		        float Sn[], 
		        float w[],
		        int   order
			);
void bw_expand_lsps(float lsp[], int order);
void locate_lsps_jnd_steps(float lsp[], int order);

#endif
