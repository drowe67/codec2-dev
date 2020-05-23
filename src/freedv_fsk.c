/*---------------------------------------------------------------------------*\

  FILE........: freedv_fsk.c
  AUTHOR......: David Rowe
  DATE CREATED: May 2020

  Functions that implement the FreeDV modes that use the FSK modem.

\*---------------------------------------------------------------------------*/

#include <assert.h>
#include <stdlib.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "fsk.h"
#include "fmfsk.h"
#include "codec2.h"
#include "freedv_vhf_framing.h"
#include "varicode.h"
#include "freedv_api.h"
#include "freedv_api_internal.h"
#include "comp_prim.h"
#include "debug_alloc.h"

void freedv_2400a_open(struct freedv *f) {
    f->n_protocol_bits = 20;
    f->deframer = fvhff_create_deframer(FREEDV_VHF_FRAME_A,0);
    assert(f->deframer != NULL);  
    f->fsk = fsk_create_hbr(48000,1200,4,10,FSK_DEFAULT_NSYM,1200,1200);
    assert(f->fsk != NULL);
        
    /* Note: fsk expects tx/rx bits as an array of uint8_ts, not ints */
    f->tx_bits = (int*)MALLOC(f->fsk->Nbits*sizeof(uint8_t));        
    assert(f->tx_bits != NULL);

    f->n_nom_modem_samples = f->fsk->N;
    f->n_max_modem_samples = f->fsk->N + (f->fsk->Ts);
    f->n_nat_modem_samples = f->fsk->N;
    f->nin = f->nin_prev = fsk_nin(f->fsk);
    f->modem_sample_rate = 48000;
    f->modem_symbol_rate = 1200;

    f->speech_sample_rate = FREEDV_FS_8000;
    f->codec2 = codec2_create(CODEC2_MODE_1300); assert(f->codec2 != NULL);
    f->n_speech_samples = codec2_samples_per_frame(f->codec2);

    f->n_codec_frames = 1;
    f->bits_per_codec_frame = codec2_bits_per_frame(f->codec2);
    f->bits_per_modem_frame = f->bits_per_codec_frame;
    int n_packed_bytes = (f->bits_per_codec_frame + 7)/8;
    f->tx_payload_bits = MALLOC(n_packed_bytes); assert(f->tx_payload_bits != NULL);  
    f->rx_payload_bits = MALLOC(n_packed_bytes); assert(f->rx_payload_bits != NULL);  
}

void freedv_2400b_open(struct freedv *f) {
    f->n_protocol_bits = 20;
    f->deframer = fvhff_create_deframer(FREEDV_VHF_FRAME_A,1);
    assert(f->deframer != NULL);
        
    f->fmfsk = fmfsk_create(48000,2400);
    assert (f->fmfsk != NULL);

    /* Note: fsk expects tx/rx bits as an array of uint8_ts, not ints */
    f->tx_bits = (int*)MALLOC(f->fmfsk->nbit*sizeof(uint8_t));
    assert(f->tx_bits != NULL);
    f->n_nom_modem_samples = f->fmfsk->N;
    f->n_max_modem_samples = f->fmfsk->N + (f->fmfsk->Ts);
    f->n_nat_modem_samples = f->fmfsk->N;
    f->nin = f->nin_prev = fmfsk_nin(f->fmfsk);
    f->modem_sample_rate = 48000;

    f->speech_sample_rate = FREEDV_FS_8000;
    f->codec2 = codec2_create(CODEC2_MODE_1300); assert(f->codec2 != NULL);
    f->n_speech_samples = codec2_samples_per_frame(f->codec2);

    f->n_codec_frames = 1;
    f->bits_per_codec_frame = codec2_bits_per_frame(f->codec2);
    f->bits_per_modem_frame = f->bits_per_codec_frame;
    int n_packed_bytes = (f->bits_per_codec_frame + 7)/8;
    f->tx_payload_bits = MALLOC(n_packed_bytes); assert(f->tx_payload_bits != NULL);  
    f->rx_payload_bits = MALLOC(n_packed_bytes); assert(f->rx_payload_bits != NULL);  
}

void freedv_800xa_open(struct freedv *f) {
    f->deframer = fvhff_create_deframer(FREEDV_HF_FRAME_B,0);
    assert(f->deframer != NULL);
    f->fsk = fsk_create_hbr(8000,400,4,10,32,800,400);
    assert(f->fsk != NULL);
        
    f->tx_bits = (int*)MALLOC(f->fsk->Nbits*sizeof(uint8_t));
    assert(f->fsk != NULL);
        
    f->n_nom_modem_samples = f->fsk->N;
    f->n_max_modem_samples = f->fsk->N + (f->fsk->Ts);
    f->n_nat_modem_samples = f->fsk->N;
    f->nin = f->nin_prev = fsk_nin(f->fsk);
    f->modem_sample_rate = 8000;
    f->modem_symbol_rate = 400;
    fsk_stats_normalise_eye(f->fsk, 0);

    f->codec2 = codec2_create(CODEC2_MODE_700C); assert(f->codec2 != NULL);
    f->speech_sample_rate = FREEDV_FS_8000;
    f->n_speech_samples = 2*codec2_samples_per_frame(f->codec2);

    f->bits_per_codec_frame = codec2_bits_per_frame(f->codec2);
    f->bits_per_modem_frame = f->n_codec_frames*f->bits_per_codec_frame;
    int n_packed_bytes = (f->bits_per_codec_frame + 7)/8;
    f->tx_payload_bits = MALLOC(n_packed_bytes); assert(f->tx_payload_bits != NULL);  
    f->rx_payload_bits = MALLOC(n_packed_bytes); assert(f->rx_payload_bits != NULL);  
}

/* TX routines for 2400 FSK modes, after codec2 encoding */
void freedv_tx_fsk_voice(struct freedv *f, short mod_out[]) {
    int  i;
    float *tx_float; /* To hold on to modulated samps from fsk/fmfsk */
    uint8_t vc_bits[2]; /* Varicode bits for 2400 framing */
    uint8_t proto_bits[3]; /* Prococol bits for 2400 framing */

    /* Frame for 2400A/B */
    if(FDV_MODE_ACTIVE( FREEDV_MODE_2400A, f->mode) || FDV_MODE_ACTIVE( FREEDV_MODE_2400B, f->mode)){
        /* Get varicode bits for TX and possibly ask for a new char */
        /* 2 bits per 2400A/B frame, so this has to be done twice */
        for(i=0;i<2;i++){
            if (f->nvaricode_bits) {
                vc_bits[i] = f->tx_varicode_bits[f->varicode_bit_index++];
                f->nvaricode_bits--;
            }

            if (f->nvaricode_bits == 0) {
                /* get new char and encode */
                char s[2];
                if (f->freedv_get_next_tx_char != NULL) {
                    s[0] = (*f->freedv_get_next_tx_char)(f->callback_state);
                    f->nvaricode_bits = varicode_encode(f->tx_varicode_bits, s, VARICODE_MAX_BITS, 1, 1);
                    f->varicode_bit_index = 0;
                }
            }
        }

        /* If the API user hasn't set up message callbacks, don't bother with varicode bits */
        if(f->freedv_get_next_proto != NULL){
            (*f->freedv_get_next_proto)(f->proto_callback_state,(char*)proto_bits);
            fvhff_frame_bits(FREEDV_VHF_FRAME_A,(uint8_t*)(f->tx_bits),f->tx_payload_bits,proto_bits,vc_bits);
        }else if(f->freedv_get_next_tx_char != NULL){
            fvhff_frame_bits(FREEDV_VHF_FRAME_A,(uint8_t*)(f->tx_bits),f->tx_payload_bits,NULL,vc_bits);
        }else {
            fvhff_frame_bits(FREEDV_VHF_FRAME_A,(uint8_t*)(f->tx_bits),f->tx_payload_bits,NULL,NULL);
        }
    /* Frame for 800XA */
    }else if(FDV_MODE_ACTIVE( FREEDV_MODE_800XA, f->mode)){
        fvhff_frame_bits(FREEDV_HF_FRAME_B,(uint8_t*)(f->tx_bits),f->tx_payload_bits,NULL,NULL);
    }

    /* Allocate floating point buffer for FSK mod */
    tx_float = alloca(sizeof(float)*f->n_nom_modem_samples);

    /* do 4fsk mod */
    if(FDV_MODE_ACTIVE( FREEDV_MODE_2400A, f->mode) || FDV_MODE_ACTIVE( FREEDV_MODE_800XA, f->mode)){
        if (f->ext_vco) {
            fsk_mod_ext_vco(f->fsk,tx_float,(uint8_t*)(f->tx_bits));
            for(i=0; i<f->n_nom_modem_samples; i++){
                mod_out[i] = (short)tx_float[i];
            }
        }
        else {
            fsk_mod(f->fsk,tx_float,(uint8_t*)(f->tx_bits));
            /* Convert float samps to short */
            for(i=0; i<f->n_nom_modem_samples; i++){
                mod_out[i] = (short)(tx_float[i]*FSK_SCALE*NORM_PWR_FSK);
            }
        }
    /* do me-fsk mod */
    }else if(FDV_MODE_ACTIVE( FREEDV_MODE_2400B, f->mode)){
        fmfsk_mod(f->fmfsk,tx_float,(uint8_t*)(f->tx_bits));
        /* Convert float samps to short */
        for(i=0; i<f->n_nom_modem_samples; i++){
            mod_out[i] = (short)(tx_float[i]*FMFSK_SCALE);
        }
    }
}

/* TX routines for 2400 FSK modes, after codec2 encoding */
void freedv_comptx_fsk_voice(struct freedv *f, COMP mod_out[]) {
    int  i;
    float *tx_float; /* To hold on to modulated samps from fsk/fmfsk */
    uint8_t vc_bits[2]; /* Varicode bits for 2400 framing */
    uint8_t proto_bits[3]; /* Prococol bits for 2400 framing */

    /* Frame for 2400A/B */
    if(FDV_MODE_ACTIVE( FREEDV_MODE_2400A, f->mode) || FDV_MODE_ACTIVE( FREEDV_MODE_2400B, f->mode)){
        /* Get varicode bits for TX and possibly ask for a new char */
        /* 2 bits per 2400A/B frame, so this has to be done twice */
        for(i=0;i<2;i++){
            if (f->nvaricode_bits) {
                vc_bits[i] = f->tx_varicode_bits[f->varicode_bit_index++];
                f->nvaricode_bits--;
            }

            if (f->nvaricode_bits == 0) {
                /* get new char and encode */
                char s[2];
                if (f->freedv_get_next_tx_char != NULL) {
                    s[0] = (*f->freedv_get_next_tx_char)(f->callback_state);
                    f->nvaricode_bits = varicode_encode(f->tx_varicode_bits, s, VARICODE_MAX_BITS, 1, 1);
                    f->varicode_bit_index = 0;
                }
            }
        }

        /* If the API user hasn't set up message callbacks, don't bother with varicode bits */
        if(f->freedv_get_next_proto != NULL){
            (*f->freedv_get_next_proto)(f->proto_callback_state,(char*)proto_bits);
            fvhff_frame_bits(FREEDV_VHF_FRAME_A,(uint8_t*)(f->tx_bits),f->tx_payload_bits,proto_bits,vc_bits);
        }else if(f->freedv_get_next_tx_char != NULL){
            fvhff_frame_bits(FREEDV_VHF_FRAME_A,(uint8_t*)(f->tx_bits),f->tx_payload_bits,NULL,vc_bits);
        }else {
            fvhff_frame_bits(FREEDV_VHF_FRAME_A,(uint8_t*)(f->tx_bits),f->tx_payload_bits,NULL,NULL);
        }
    /* Frame for 800XA */
    }else if(FDV_MODE_ACTIVE( FREEDV_MODE_800XA, f->mode)){
        fvhff_frame_bits(FREEDV_HF_FRAME_B,(uint8_t*)(f->tx_bits),f->tx_payload_bits,NULL,NULL);
    }

    /* Allocate floating point buffer for FSK mod */
    tx_float = alloca(sizeof(float)*f->n_nom_modem_samples);

    /* do 4fsk mod */
    if(FDV_MODE_ACTIVE( FREEDV_MODE_2400A, f->mode) || FDV_MODE_ACTIVE( FREEDV_MODE_800XA, f->mode)){
        fsk_mod_c(f->fsk,mod_out,(uint8_t*)(f->tx_bits));
        /* Convert float samps to short */
        for(i=0; i<f->n_nom_modem_samples; i++){
        	mod_out[i] = fcmult(NORM_PWR_FSK,mod_out[i]);
        }
    /* do me-fsk mod */
    }else if(FDV_MODE_ACTIVE( FREEDV_MODE_2400B, f->mode)){
        fmfsk_mod(f->fmfsk,tx_float,(uint8_t*)(f->tx_bits));
        /* Convert float samps to short */
        for(i=0; i<f->n_nom_modem_samples; i++){
            mod_out[i].real = (tx_float[i]);
        }
    }
}

/* TX routines for 2400 FSK modes, data channel */
void freedv_tx_fsk_data(struct freedv *f, short mod_out[]) {
    int  i;
    float *tx_float; /* To hold on to modulated samps from fsk/fmfsk */
    
    if (FDV_MODE_ACTIVE( FREEDV_MODE_800XA, f->mode))
	fvhff_frame_data_bits(f->deframer, FREEDV_HF_FRAME_B,(uint8_t*)(f->tx_bits));
    else
    	fvhff_frame_data_bits(f->deframer, FREEDV_VHF_FRAME_A,(uint8_t*)(f->tx_bits));
        
    /* Allocate floating point buffer for FSK mod */
    tx_float = alloca(sizeof(float)*f->n_nom_modem_samples);
        
    /* do 4fsk mod */
    if(FDV_MODE_ACTIVE( FREEDV_MODE_2400A, f->mode) || FDV_MODE_ACTIVE( FREEDV_MODE_800XA, f->mode)){
        fsk_mod(f->fsk,tx_float,(uint8_t*)(f->tx_bits));
        /* Convert float samps to short */
        for(i=0; i<f->n_nom_modem_samples; i++){
            mod_out[i] = (short)(tx_float[i]*FSK_SCALE);
        }
    /* do me-fsk mod */
    }else if(FDV_MODE_ACTIVE( FREEDV_MODE_2400B, f->mode)){
        fmfsk_mod(f->fmfsk,tx_float,(uint8_t*)(f->tx_bits));
        /* Convert float samps to short */
        for(i=0; i<f->n_nom_modem_samples; i++){
            mod_out[i] = (short)(tx_float[i]*FMFSK_SCALE);
        }
    }
}


// float input samples version
int freedv_comprx_fsk(struct freedv *f, COMP demod_in[]) {
    /* Varicode and protocol bits */
    uint8_t vc_bits[2];
    uint8_t proto_bits[3];
    short vc_bit;
    int i;
    int n_ascii;
    char ascii_out;
    int rx_status = 0;
    
    if(FDV_MODE_ACTIVE( FREEDV_MODE_2400A, f->mode) || FDV_MODE_ACTIVE( FREEDV_MODE_800XA, f->mode)){        
	fsk_demod(f->fsk,(uint8_t*)f->tx_bits,demod_in);
        f->nin = fsk_nin(f->fsk);
        float EbNodB = f->fsk->stats->snr_est;           /* fsk demod actually estimates Eb/No     */
        f->snr_est = EbNodB + 10.0*log10f(800.0/3000.0); /* so convert to SNR Rb=800, noise B=3000 */
    } else{      
        /* 2400B needs real input samples */
        int n = fmfsk_nin(f->fmfsk);
        float demod_in_float[n];
        for(i=0; i<n; i++) {
            demod_in_float[i] = demod_in[i].real;
        }
        fmfsk_demod(f->fmfsk,(uint8_t*)f->tx_bits,demod_in_float);
        f->nin = fmfsk_nin(f->fmfsk);
    }
    
    if(fvhff_deframe_bits(f->deframer,f->rx_payload_bits,proto_bits,vc_bits,(uint8_t*)f->tx_bits)){
        /* Decode varicode text */
        for(i=0; i<2; i++){
            /* Note: deframe_bits spits out bits in uint8_ts while varicode_decode expects shorts */
            vc_bit = vc_bits[i];
            n_ascii = varicode_decode(&f->varicode_dec_states, &ascii_out, &vc_bit, 1, 1);
            if (n_ascii && (f->freedv_put_next_rx_char != NULL)) {
                (*f->freedv_put_next_rx_char)(f->callback_state, ascii_out);
            }
        }
        /* Pass proto bits on down if callback is present */
        if( f->freedv_put_next_proto != NULL){
            (*f->freedv_put_next_proto)(f->proto_callback_state,(char*)proto_bits);
        }

        rx_status = RX_SYNC | RX_BITS;
    } 
    f->sync = f->deframer->state;
    f->stats.sync = f->deframer->state;

    return rx_status;
}


int freedv_floatrx(struct freedv *f, short speech_out[], float demod_in[]) {
    assert(f != NULL);
    int  i;
    int nin = freedv_nin(f);    
    
    assert(nin <= f->n_max_modem_samples);
    
    COMP rx_fdm[f->n_max_modem_samples];
    for(i=0; i<nin; i++) {
        rx_fdm[i].real = demod_in[i];
        rx_fdm[i].imag = 0;
    }

    return freedv_comprx(f, speech_out, rx_fdm);
}
