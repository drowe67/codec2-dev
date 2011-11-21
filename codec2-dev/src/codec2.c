/*---------------------------------------------------------------------------*\

  FILE........: codec2.c
  AUTHOR......: David Rowe
  DATE CREATED: 21/8/2010

  Codec2 fully quantised encoder and decoder functions.  If you want use 
  codec2, the codec2_xxx functions are for you.

\*---------------------------------------------------------------------------*/

/*
  Copyright (C) 2010 David Rowe

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

#include "defines.h"
#include "sine.h"
#include "nlp.h"
#include "dump.h"
#include "lpc.h"
#include "quantise.h"
#include "phase.h"
#include "interp.h"
#include "postfilter.h"
#include "codec2.h"
#include "lsp.h"

struct CODEC2 {
    int    mode;
    float  w[M];	        /* time domain hamming window                */
    COMP   W[FFT_ENC];	        /* DFT of w[]                                */
    float  Pn[2*N];	        /* trapezoidal synthesis window              */
    float  Sn[M];               /* input speech                              */
    float  hpf_states[2];       /* high pass filter states                   */
    void  *nlp;                 /* pitch predictor states                    */
    float  Sn_[2*N];	        /* synthesised output speech                 */
    float  ex_phase;            /* excitation model phase track              */
    float  bg_est;              /* background noise estimate for post filter */
    float  prev_Wo;             /* previous frame's pitch estimate           */
    MODEL  prev_model;          /* previous frame's model parameters         */
    float  prev_lsps_[LPC_ORD]; /* previous frame's LSPs                     */
    float  prev_energy;         /* previous frame's LPC energy               */
};

/*---------------------------------------------------------------------------*\
                                                       
                             FUNCTION HEADERS

\*---------------------------------------------------------------------------*/

void analyse_one_frame(struct CODEC2 *c2, MODEL *model, short speech[]);
void synthesise_one_frame(struct CODEC2 *c2, short speech[], MODEL *model,
			  float ak[]);
void codec2_encode_2500(struct CODEC2 *c2, unsigned char * bits, short speech[]);
void codec2_decode_2500(struct CODEC2 *c2, short speech[], const unsigned char * bits);
void codec2_encode_1400(struct CODEC2 *c2, unsigned char * bits, short speech[]);
void codec2_decode_1400(struct CODEC2 *c2, short speech[], const unsigned char * bits);

/*---------------------------------------------------------------------------*\
                                                       
                                FUNCTIONS

\*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*\
                                                       
  FUNCTION....: codec2_create	     
  AUTHOR......: David Rowe			      
  DATE CREATED: 21/8/2010 

  Create and initialise an instance of the codec.  Returns a pointer
  to the codec states or NULL on failure.  One set of states is
  sufficient for a full duuplex codec (i.e. an encoder and decoder).
  You don't need separate states for encoders and decoders.  See
  c2enc.c and c2dec.c for examples.

\*---------------------------------------------------------------------------*/

struct CODEC2 *codec2_create(int mode)
{
    struct CODEC2 *c2;
    int            i,l;

    c2 = (struct CODEC2*)malloc(sizeof(struct CODEC2));
    if (c2 == NULL)
	return NULL;
    
    assert((mode == CODEC2_MODE_2500) || (mode == CODEC2_MODE_1400));
    c2->mode = mode;
    for(i=0; i<M; i++)
	c2->Sn[i] = 1.0;
    c2->hpf_states[0] = c2->hpf_states[1] = 0.0;
    for(i=0; i<2*N; i++)
	c2->Sn_[i] = 0;
    make_analysis_window(c2->w,c2->W);
    make_synthesis_window(c2->Pn);
    quantise_init();
    c2->prev_Wo = 0.0;
    c2->bg_est = 0.0;
    c2->ex_phase = 0.0;

    for(l=1; l<MAX_AMP; l++)
	c2->prev_model.A[l] = 0.0;
    c2->prev_model.Wo = TWO_PI/P_MAX;
    c2->prev_model.L = PI/c2->prev_model.Wo;
    c2->prev_model.voiced = 0;

    for(i=0; i<LPC_ORD; i++) {
      c2->prev_lsps_[i] = i*PI/(LPC_ORD+1);
    }
    c2->prev_energy = 1;

    c2->nlp = nlp_create();
    if (c2->nlp == NULL) {
	free (c2);
	return NULL;
    }

    return c2;
}

/*---------------------------------------------------------------------------*\
                                                       
  FUNCTION....: codec2_destroy	     
  AUTHOR......: David Rowe			      
  DATE CREATED: 21/8/2010 

  Destroy an instance of the codec.

\*---------------------------------------------------------------------------*/

void codec2_destroy(struct CODEC2 *c2)
{
    assert(c2 != NULL);
    nlp_destroy(c2->nlp);
    free(c2);
}

/*---------------------------------------------------------------------------*\
                                                       
  FUNCTION....: codec2_bits_per_frame     
  AUTHOR......: David Rowe			      
  DATE CREATED: Nov 14 2011

  Returns the number of bits per frame.

\*---------------------------------------------------------------------------*/

int codec2_bits_per_frame(struct CODEC2 *c2) {
    if (c2->mode == CODEC2_MODE_2500)
	return 50;
    if  (c2->mode == CODEC2_MODE_1400)
	return 56;

    return 0; /* shouldn't get here */
}


/*---------------------------------------------------------------------------*\
                                                       
  FUNCTION....: codec2_samples_per_frame     
  AUTHOR......: David Rowe			      
  DATE CREATED: Nov 14 2011

  Returns the number of bits per frame.

\*---------------------------------------------------------------------------*/

int codec2_samples_per_frame(struct CODEC2 *c2) {
    if (c2->mode == CODEC2_MODE_2500)
	return 160;
    if  (c2->mode == CODEC2_MODE_1400)
	return 320;

    return 0; /* shouldnt get here */
}

void codec2_encode(struct CODEC2 *c2, unsigned char *bits, short speech[])
{
    assert(c2 != NULL);
    assert((c2->mode == CODEC2_MODE_2500) || (c2->mode == CODEC2_MODE_1400));

    if (c2->mode == CODEC2_MODE_2500)
	codec2_encode_2500(c2, bits, speech);
    else
	codec2_encode_1400(c2, bits, speech);
}

void codec2_decode(struct CODEC2 *c2, short speech[], const unsigned char *bits)
{
    assert(c2 != NULL);
    assert((c2->mode == CODEC2_MODE_2500) || (c2->mode == CODEC2_MODE_1400));

    if (c2->mode == CODEC2_MODE_2500)
	codec2_decode_2500(c2, speech, bits);
    else
	codec2_decode_1400(c2, speech, bits);
}

/*---------------------------------------------------------------------------*\
                                                       
  FUNCTION....: codec2_encode_2500	     
  AUTHOR......: David Rowe			      
  DATE CREATED: 21/8/2010 

  Encodes 160 speech samples (20ms of speech) into 50 bits.  

  The codec2 algorithm actually operates internally on 10ms (80
  sample) frames, so we run the encoding algorithm twice.  On the
  first frame we just send the voicing bit.  One the second frame we
  send all model parameters.

  The bit allocation is:

    Parameter                      bits/frame
    --------------------------------------
    Harmonic magnitudes (LSPs)     36
    Energy                          5
    Wo (fundamental frequnecy)      7
    Voicing (10ms update)           2
    TOTAL                          50
 
\*---------------------------------------------------------------------------*/

void codec2_encode_2500(struct CODEC2 *c2, unsigned char * bits, short speech[])
{
    MODEL   model;
    int     voiced1, voiced2;
    float   lsps[LPC_ORD];
    float   ak[LPC_ORD+1];
    float   e;
    int     lsp_indexes[LPC_ORD];
    int     energy_index;
    int     Wo_index;
    int     i;
    unsigned int nbit = 0;

    assert(c2 != NULL);

    memset(bits, '\0', ((codec2_bits_per_frame(c2) + 7) / 8));

    /* first 10ms analysis frame - we just want voicing */

    analyse_one_frame(c2, &model, speech);
    voiced1 = model.voiced;

    /* second 10ms analysis frame */

    analyse_one_frame(c2, &model, &speech[N]);
    voiced2 = model.voiced;
    
    Wo_index = encode_Wo(model.Wo);

    e = speech_to_uq_lsps(lsps, ak, c2->Sn, c2->w, LPC_ORD);
    encode_lsps_scalar(lsp_indexes, lsps, LPC_ORD);
    energy_index = encode_energy(e);
    //for(i=0; i<LPC_ORD; i++)
    //	printf("lsp_indexes: %d lsps: %2.3f\n", lsp_indexes[i], lsps[i]);
    //exit(0);

    pack(bits, &nbit, Wo_index, WO_BITS);
    for(i=0; i<LSP_SCALAR_INDEXES; i++) {
	pack(bits, &nbit, lsp_indexes[i], lsp_bits(i));
    }
    pack(bits, &nbit, energy_index, E_BITS);
    pack(bits, &nbit, voiced1, 1);
    pack(bits, &nbit, voiced2, 1);
    printf("v2: %d  v1: %d\n", voiced2, voiced1);
    assert(nbit == codec2_bits_per_frame(c2));
}


/*---------------------------------------------------------------------------*\
                                                       
  FUNCTION....: codec2_decode_2500	     
  AUTHOR......: David Rowe			      
  DATE CREATED: 21/8/2010 

  Decodes frames of 50 bits into 160 samples (20ms) of speech.

\*---------------------------------------------------------------------------*/

void codec2_decode_2500(struct CODEC2 *c2, short speech[], const unsigned char * bits)
{
    MODEL   model;
    int     voiced1, voiced2;
    int     lsp_indexes[LPC_ORD];
    float   lsps_[LPC_ORD];
    int     energy_index;
    float   energy;
    float   snr;
    int     Wo_index;
    float   ak[LPC_ORD+1];
    float   ak_interp[LPC_ORD+1];
    int     i;
    unsigned int nbit = 0;
    MODEL   model_interp;
    static  int frames;

    printf("frame: %d\n", frames+=2);
    assert(c2 != NULL);
    
    /* unpack bit stream to integer codes */

    Wo_index = unpack(bits, &nbit, WO_BITS);
    for(i=0; i<LSP_SCALAR_INDEXES; i++) {
	lsp_indexes[i] = unpack(bits, &nbit, lsp_bits(i));
    }
    energy_index = unpack(bits, &nbit, E_BITS);
    voiced1 = unpack(bits, &nbit, 1);
    voiced2 = unpack(bits, &nbit, 1);
    assert(nbit == codec2_bits_per_frame(c2));

    /* decode integer codes to model parameters */

    model.Wo = decode_Wo(Wo_index);
    model.L = PI/model.Wo;
    memset(&model.A, 0, (model.L+1)*sizeof(model.A[0]));

    /* decode even frame LSPs and model amplitudes */

    decode_lsps_scalar(lsps_, lsp_indexes, LPC_ORD);
    bw_expand_lsps(lsps_, LPC_ORD);
    lsp_to_lpc(lsps_, ak, LPC_ORD);
    energy = decode_energy(energy_index);
    aks_to_M2(ak, LPC_ORD, &model, energy, &snr, 1); 
    apply_lpc_correction(&model);

    printf("Wo: %1.5f  L: %d e: %3.2f v2: %d\n", 
	   model.Wo, model.L, energy, voiced2 );
    //for(i=0; i<LPC_ORD; i++)
    //	printf("lsp_indexes: %d lsp_: %2.3f prev_lsp_: %2.3f\n", 
    //	       lsp_indexes[i], lsps_[i], c2->prev_lsps_[i]);
    //printf("ak: ");
    //for(i=0; i<LPC_ORD; i++)
    //	printf("%2.3f  ", ak[i]);
    printf("Am: ");
    for(i=0; i<5; i++)
    	printf("%2.3f  ", model.A[i]);
    printf("\n");
    
    /* interpolate odd frame model parameters from adjacent frames */

    model.voiced = voiced2;
    model_interp.voiced = voiced1;
    model_interp.Wo = P_MAX/2;
    memset(&model_interp.A, 0, MAX_AMP*sizeof(model_interp.A[0]));

    interpolate_lsp(&model_interp, &c2->prev_model, &model,
    		    c2->prev_lsps_, c2->prev_energy, lsps_, energy, ak_interp);
    apply_lpc_correction(&model_interp);
    printf("Wo: %1.5f  L: %d prev_e: %3.2f v1: %d pv: %d\n", 
	   model_interp.Wo, model_interp.L, c2->prev_energy, voiced1,
	   c2->prev_model.voiced);
    //printf("ak_interp: ");
    //for(i=0; i<LPC_ORD; i++)
    //	printf("%2.3f  ", ak_interp[i]);
    //printf("\n");
    printf("Am: ");
    for(i=0; i<5; i++)
    	printf("%2.3f  ", model_interp.A[i]);
    printf("\n");
    //if (frames == 6)
    //	exit(0);

    /* synthesise two 10ms frames */

    synthesise_one_frame(c2, speech, &model_interp, ak_interp);
    synthesise_one_frame(c2, &speech[N], &model, ak);

    /* update memories (decode states) for next time */

    memcpy(&c2->prev_model, &model, sizeof(MODEL));
    memcpy(c2->prev_lsps_, lsps_, sizeof(lsps_));
    c2->prev_energy = energy;
}


/*---------------------------------------------------------------------------*\
                                                       
  FUNCTION....: codec2_encode_1400	     
  AUTHOR......: David Rowe			      
  DATE CREATED: Nov 14 2011 

  Encodes 320 speech samples (50ms of speech) into 56 bits.  

  The codec2 algorithm actually operates internally on 10ms (80
  sample) frames, so we run the encoding algorithm for times:

  frame 0: just send voicing bit
  frame 1: full quantisation of LSPs and Wo
  frame 2: just send voicing bit
  frame 3: delta-time quantisation of LSPs and Wo

  The bit allocation is:

    Parameter                      frame 2  frame 4   Total
    -------------------------------------------------------
    Harmonic magnitudes (LSPs)     25        7        32
    Energy                          5        5        10
    Wo (fundamental frequnecy)      7        3        10
    Voicing (10ms update)           2        2         4
    TOTAL                          39       17        56
 
\*---------------------------------------------------------------------------*/

void codec2_encode_1400(struct CODEC2 *c2, unsigned char * bits, short speech[])
{
    MODEL   model;
    float   lsps[LPC_ORD], prev_lsps_[LPC_ORD];
    float   ak[LPC_ORD+1];
    float   e;
    int     voiced1, voiced2, voiced3, voiced4;
    int     lsp_indexes[LPC_ORD];
    int     energy_index;
    int     Wo_index, delta_Wo_index;
    int     i;
    unsigned int nbit = 0;
    float   prev_Wo;

    assert(c2 != NULL);

    memset(bits, '\0', codec2_bits_per_frame(c2)/8);

    /* frame 1: - we just want voicing -------------------------------- */

    analyse_one_frame(c2, &model, speech);
    voiced1 = model.voiced;

    /* frame 2: - full LSP and Wo ------------------------------------- */

    analyse_one_frame(c2, &model, &speech[N]);
    voiced2 = model.voiced;
    
    Wo_index = encode_Wo(model.Wo);

    e = speech_to_uq_lsps(lsps, ak, c2->Sn, c2->w, LPC_ORD);
    encode_lsps_diff_freq_vq(lsp_indexes, lsps, LPC_ORD);
    energy_index = encode_energy(e);

    pack(bits, &nbit, Wo_index, WO_BITS);
    for(i=0; i<LSP_DIFF_FREQ_INDEXES; i++) {
	pack(bits, &nbit, lsp_indexes[i], lspd_bits(i));
    }
    pack(bits, &nbit, energy_index, E_BITS);
    pack(bits, &nbit, voiced1, 1);
    pack(bits, &nbit, voiced2, 1);

    /* decode LSPs and Wo for use in frame 4 */

    decode_lsps_diff_freq_vq(prev_lsps_, lsp_indexes, LPC_ORD);
    prev_Wo = decode_Wo(Wo_index);

    /* frame 3: - we just want voicing --------------------------------- */

    analyse_one_frame(c2, &model, &speech[2*N]);
    voiced3 = model.voiced;

    /* frame 4: - delta-time and delta Wo -----------------------------  */

    analyse_one_frame(c2, &model, &speech[3*N]);
    voiced4 = model.voiced;
    
    delta_Wo_index =  encode_Wo_dt(model.Wo, prev_Wo);
  
    e = speech_to_uq_lsps(lsps, ak, c2->Sn, c2->w, LPC_ORD);
    encode_lsps_diff_time_vq(lsp_indexes, lsps, prev_lsps_, LPC_ORD);
    energy_index = encode_energy(e);

    pack(bits, &nbit, delta_Wo_index, WO_DT_BITS);
    pack(bits, &nbit, *lsp_indexes, LSP_DIFF_TIME_BITS);
    pack(bits, &nbit, energy_index, E_BITS);
    pack(bits, &nbit, voiced3, 1);
    pack(bits, &nbit, voiced4, 1);
    
    assert(nbit == codec2_bits_per_frame(c2));
}


/*---------------------------------------------------------------------------*\
                                                       
  FUNCTION....: codec2_decode_1400	     
  AUTHOR......: David Rowe			      
  DATE CREATED: 16 Nov 2011

  Decodes frames of 56 bits into 320 samples (40ms) of speech.

\*---------------------------------------------------------------------------*/

void codec2_decode_1400(struct CODEC2 *c2, short speech[], const unsigned char * bits)
{
    MODEL   model;
    int     voiced1, voiced2, voiced3, voiced4;
    int     lsp_indexes[LPC_ORD];
    float   lsps_[LPC_ORD];
    int     energy_index;
    float   energy;
    float   snr;
    int     Wo_index, delta_Wo_index;
    float   ak[LPC_ORD+1];
    float   ak_interp[LPC_ORD+1];
    int     i;
    unsigned int nbit = 0;
    MODEL   model_interp;

    assert(c2 != NULL);

    /* unpack frame 1 & 2 bit stream to integer codes */

    Wo_index = unpack(bits, &nbit, WO_BITS);
    for(i=0; i<LSP_DIFF_FREQ_INDEXES; i++) {
	lsp_indexes[i] = unpack(bits, &nbit, lspd_bits(i));
    }
    energy_index = unpack(bits, &nbit, E_BITS);
    voiced1 = unpack(bits, &nbit, 1);
    voiced2 = unpack(bits, &nbit, 1);
 
    /* decode integer codes to model parameters */

    model.Wo = decode_Wo(Wo_index);
    model.L = PI/model.Wo;
    memset(&model.A, 0, (model.L+1)*sizeof(model.A[0]));

    /* decode frame 2 LSPs and model amplitudes */

    decode_lsps_diff_freq_vq(lsps_, lsp_indexes, LPC_ORD);
    bw_expand_lsps(lsps_, LPC_ORD);
    lsp_to_lpc(lsps_, ak, LPC_ORD);
    energy = decode_energy(energy_index);
    aks_to_M2(ak, LPC_ORD, &model, energy, &snr, 1); 
    apply_lpc_correction(&model);

    /* interpolate frame 1 model parameters from adjacent frames */

    model.voiced = voiced2;
    model_interp.voiced = voiced1;
    model_interp.Wo = P_MAX/2;
    memset(&model_interp.A, 0, MAX_AMP*sizeof(model_interp.A[0]));

    interpolate_lsp(&model_interp, &c2->prev_model, &model,
    		    c2->prev_lsps_, c2->prev_energy, lsps_, energy, ak_interp);
    apply_lpc_correction(&model_interp);

    /* synthesise frame 1 and frame 2 10ms frames */

    synthesise_one_frame(c2, speech, &model_interp, ak_interp);
    synthesise_one_frame(c2, &speech[N], &model, ak);

    /* update memories (decode states) for next time */

    memcpy(&c2->prev_model, &model, sizeof(MODEL));
    memcpy(c2->prev_lsps_, lsps_, sizeof(lsps_));
    c2->prev_energy = energy;

    /*--------------------------------------------------------------------*/

    /* unpack frame 3 & 4 bit stream to integer codes */

    delta_Wo_index = unpack(bits, &nbit, WO_DT_BITS);
    *lsp_indexes = unpack(bits, &nbit, LSP_DIFF_TIME_BITS);
    energy_index = unpack(bits, &nbit, E_BITS);
    voiced3 = unpack(bits, &nbit, 1);
    voiced4 = unpack(bits, &nbit, 1);
    assert(nbit == codec2_bits_per_frame(c2));

    /* decode frame 4 LSPs and model amplitudes */

    decode_lsps_diff_time_vq(lsps_, lsp_indexes, c2->prev_lsps_, LPC_ORD);
    bw_expand_lsps(lsps_, LPC_ORD);
    lsp_to_lpc(lsps_, ak, LPC_ORD);
    energy = decode_energy(energy_index);
    aks_to_M2(ak, LPC_ORD, &model, energy, &snr, 1); 
    apply_lpc_correction(&model);

    /* interpolate frame 3 model parameters from adjacent frames */

    model.voiced = voiced4;
    model_interp.voiced = voiced3;
    model_interp.Wo = P_MAX/2;
    memset(&model_interp.A, 0, MAX_AMP*sizeof(model_interp.A[0]));

    interpolate_lsp(&model_interp, &c2->prev_model, &model,
    		    c2->prev_lsps_, c2->prev_energy, lsps_, energy, ak_interp);
    apply_lpc_correction(&model_interp);

    /* synthesise frame 3 and frame 4 10ms frames */

    synthesise_one_frame(c2, &speech[2*N], &model_interp, ak_interp);
    synthesise_one_frame(c2, &speech[3*N], &model, ak);

    /* update memories (decode states) for next time */

    memcpy(&c2->prev_model, &model, sizeof(MODEL));
    memcpy(c2->prev_lsps_, lsps_, sizeof(lsps_));
    c2->prev_energy = energy;

}


/*---------------------------------------------------------------------------*\
                                                       
  FUNCTION....: synthesise_one_frame()	     
  AUTHOR......: David Rowe			      
  DATE CREATED: 23/8/2010 

  Synthesise 80 speech samples (10ms) from model parameters.

\*---------------------------------------------------------------------------*/

void synthesise_one_frame(struct CODEC2 *c2, short speech[], MODEL *model, float ak[])
{
    int     i;

    phase_synth_zero_order(model, ak, &c2->ex_phase, LPC_ORD);
    postfilter(model, &c2->bg_est);
    synthesise(c2->Sn_, model, c2->Pn, 1);

    for(i=0; i<N; i++) {
	if (c2->Sn_[i] > 32767.0)
	    speech[i] = 32767;
	else if (c2->Sn_[i] < -32767.0)
	    speech[i] = -32767;
	else
	    speech[i] = c2->Sn_[i];
    }

}

/*---------------------------------------------------------------------------*\
                                                       
  FUNCTION....: analyse_one_frame()   
  AUTHOR......: David Rowe			      
  DATE CREATED: 23/8/2010 

  Extract sinusoidal model parameters from 80 speech samples (10ms of
  speech).
 
\*---------------------------------------------------------------------------*/

void analyse_one_frame(struct CODEC2 *c2, MODEL *model, short speech[])
{
    COMP    Sw[FFT_ENC];
    COMP    Sw_[FFT_ENC];
    COMP    Ew[FFT_ENC];
    float   pitch, snr;
    int     i;

    /* Read input speech */

    for(i=0; i<M-N; i++)
      c2->Sn[i] = c2->Sn[i+N];
    for(i=0; i<N; i++)
      c2->Sn[i+M-N] = speech[i];

    dft_speech(Sw, c2->Sn, c2->w);

    /* Estimate pitch */

    nlp(c2->nlp,c2->Sn,N,M,P_MIN,P_MAX,&pitch,Sw, &c2->prev_Wo);
    model->Wo = TWO_PI/pitch;
    model->L = PI/model->Wo;

    /* estimate model parameters */

    two_stage_pitch_refinement(model, Sw);
    estimate_amplitudes(model, Sw, c2->W);
    snr = est_voicing_mbe(model, Sw, c2->W, Sw_, Ew, c2->prev_Wo);
    printf("snr %3.2f  v: %d  Wo: %f prev_Wo: %f\n", 
	   snr, model->voiced, model->Wo, c2->prev_Wo);
    c2->prev_Wo = model->Wo;
}
