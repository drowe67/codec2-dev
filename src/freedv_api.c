/*---------------------------------------------------------------------------*\

  FILE........: freedv_api.c
  AUTHOR......: David Rowe
  DATE CREATED: August 2014

  Library of API functions that implement FreeDV "modes", useful for
  embedding FreeDV in other programs.

\*---------------------------------------------------------------------------*/

/*
  Copyright (C) 2014 David Rowe

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
#include <stdlib.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#ifdef TT
#if defined(__APPLE__)
#include <malloc/malloc.h>
#elif defined(__OpenBSD__) || defined(__NetBSD__) || defined(__FreeBSD__)
#include <sys/malloc.h>
#else
#include <malloc.h>
#endif /* __APPLE__ */
#endif

#include "fsk.h"
#include "fmfsk.h"
#include "codec2.h"
#include "codec2_fdmdv.h"
#include "fdmdv_internal.h"
#include "golay23.h"
#include "varicode.h"
#include "freedv_api.h"
#include "freedv_api_internal.h"
#include "freedv_vhf_framing.h"
#include "comp_prim.h"
#include "filter.h"

#include "codec2_ofdm.h"
#include "ofdm_internal.h"
#include "mpdecode_core.h"
#include "gp_interleaver.h"
#include "interldpc.h"

#include "debug_alloc.h"

#define VERSION     12    /* The API version number.  The first version
                           is 10.  Increment if the API changes in a
                           way that would require changes by the API
                           user. */
/*
 * Version 10   Initial version August 2, 2015.
 * Version 11   September 2015
 *              Added: freedv_zero_total_bit_errors(), freedv_get_sync()
 *              Changed all input and output sample rates to 8000 sps.  Rates for FREEDV_MODE_700 and 700B were 7500.
 * Version 12   August 2018
 *              Added OFDM configuration switch structure
 */

/* experimentally derived fudge factors to normalise power across modes */

#define NORM_PWR_COHPSK  1.74   
#define NORM_PWR_FSK     0.193 
#define NORM_PWR_OFDM    1.00

static struct OFDM_CONFIG *ofdm_config;

static int ofdm_bitsperframe;
static int ofdm_nuwbits;
static int ofdm_ntxtbits;

static char *statemode[] = {
    "search",
    "trial",
    "synced"
};

/*---------------------------------------------------------------------------*\

  FUNCTION....: freedv_open
  AUTHOR......: David Rowe
  DATE CREATED: 3 August 2014

  Call this first to initialise.  Returns NULL if initialisation fails
  (e.g. out of memory or mode not supported).

\*---------------------------------------------------------------------------*/

struct freedv *freedv_open(int mode) {
    return freedv_open_advanced(mode, NULL);
}

struct freedv *freedv_open_advanced(int mode, struct freedv_advanced *adv) {
    struct freedv *f;
    int            Nc, codec2_mode, nbit, nbyte = 0;
    
    if (false == (FDV_MODE_ACTIVE( FREEDV_MODE_1600,mode) || FDV_MODE_ACTIVE( FREEDV_MODE_700,mode) || 
		FDV_MODE_ACTIVE( FREEDV_MODE_700B,mode) || FDV_MODE_ACTIVE( FREEDV_MODE_2400A,mode) || 
		FDV_MODE_ACTIVE( FREEDV_MODE_2400B,mode) || FDV_MODE_ACTIVE( FREEDV_MODE_800XA,mode) || 
		FDV_MODE_ACTIVE( FREEDV_MODE_700C,mode) || FDV_MODE_ACTIVE( FREEDV_MODE_700D,mode)  ||
                  FDV_MODE_ACTIVE( FREEDV_MODE_2020,mode)) )
    {
        return NULL;
    }
    
    f = (struct freedv*)MALLOC(sizeof(struct freedv));
    if (f == NULL)
        return NULL;

    f->mode = mode;
    f->verbose = 0;
    f->test_frames = f->smooth_symbols = 0;
    f->freedv_put_error_pattern = NULL;
    f->error_pattern_callback_state = NULL;
    f->n_protocol_bits = 0;
    f->frames = 0;
    f->speech_sample_rate = FS_VOICE_8K;
    f->stats.snr_est = 0.0;
    f->tx_bits = f->rx_bits = NULL;
    
    /* -----------------------------------------------------------------------------------------------*\
    |                     Init states for this mode, and set up samples in/out                         |
    \*------------------------------------------------------------------------------------------------*/
    
    if (FDV_MODE_ACTIVE( FREEDV_MODE_1600, f->mode)) {
        f->snr_squelch_thresh = 2.0;
        f->squelch_en = 1;
        Nc = 16;
        f->tx_sync_bit = 0;
        codec2_mode = CODEC2_MODE_1300;
        f->fdmdv = fdmdv_create(Nc);
        if (f->fdmdv == NULL)
            return NULL;
        golay23_init();
        f->nin = FDMDV_NOM_SAMPLES_PER_FRAME;
        f->n_nom_modem_samples = 2*FDMDV_NOM_SAMPLES_PER_FRAME;
        f->n_nat_modem_samples = f->n_nom_modem_samples;
        f->n_max_modem_samples = FDMDV_NOM_SAMPLES_PER_FRAME+FDMDV_MAX_SAMPLES_PER_FRAME;
        f->modem_sample_rate = FS;
        nbit = fdmdv_bits_per_frame(f->fdmdv);
        f->fdmdv_bits = (int*)MALLOC(nbit*sizeof(int));
        if (f->fdmdv_bits == NULL)
            return NULL;
        nbit = 2*fdmdv_bits_per_frame(f->fdmdv);
        f->tx_bits = (int*)CALLOC(1, nbit*sizeof(int));
        f->rx_bits = (int*)CALLOC(1, nbit*sizeof(int));
        if ((f->tx_bits == NULL) || (f->rx_bits == NULL)) {
            if (f->tx_bits != NULL) {
              FREE(f->tx_bits);
              f->tx_bits = NULL;
            }
            if (f->rx_bits != NULL) {
              FREE(f->rx_bits);
              f->rx_bits = NULL;
            }
            return NULL;
        }
        f->evenframe = 0;
        f->sz_error_pattern = fdmdv_error_pattern_size(f->fdmdv);
    }

    if (FDV_MODE_ACTIVE( FREEDV_MODE_700, mode) || FDV_MODE_ACTIVE( FREEDV_MODE_700B, mode) || FDV_MODE_ACTIVE( FREEDV_MODE_700C, mode)) {
        f->snr_squelch_thresh = 0.0;
        f->squelch_en = 0;
        switch(mode) {
        case FREEDV_MODE_700:
            codec2_mode = CODEC2_MODE_700;
            break;
        case FREEDV_MODE_700B:
            codec2_mode = CODEC2_MODE_700B;
            break;
        case FREEDV_MODE_700C:
            codec2_mode = CODEC2_MODE_700C;
            break;
        default:
            assert(0);
        }

        f->cohpsk = cohpsk_create();
        f->nin = COHPSK_NOM_SAMPLES_PER_FRAME;
        f->n_nat_modem_samples = COHPSK_NOM_SAMPLES_PER_FRAME;             // native modem samples as used by the modem
        f->n_nom_modem_samples = f->n_nat_modem_samples * FS / COHPSK_FS;  // number of samples after native samples are interpolated to 8000 sps
        f->n_max_modem_samples = COHPSK_MAX_SAMPLES_PER_FRAME * FS / COHPSK_FS + 1;
        f->modem_sample_rate = FS;                                         /* note wierd sample rate tamed by interpolator */
        f->clip = 1;
        nbit = COHPSK_BITS_PER_FRAME;
        f->tx_bits = (int*)MALLOC(nbit*sizeof(int));
        if (f->tx_bits == NULL)
            return NULL;
        f->sz_error_pattern = cohpsk_error_pattern_size();
    }
   
    if (FDV_MODE_ACTIVE( FREEDV_MODE_700D, mode) ) {
        f->snr_squelch_thresh = 0.0;
        f->squelch_en = 0;
        codec2_mode = CODEC2_MODE_700C;

        if ((ofdm_config = (struct OFDM_CONFIG *) CALLOC(1, sizeof (struct OFDM_CONFIG))) == NULL) {
            if (f->tx_bits != NULL) {
              FREE(f->tx_bits);
              f->tx_bits = NULL;
            }
            if (f->rx_bits != NULL) {
              FREE(f->rx_bits);
              f->rx_bits = NULL;
            }
            return NULL;
        }

        f->ofdm = ofdm_create(ofdm_config);
	assert(f->ofdm != NULL);
        FREE(ofdm_config);

        /* Get a copy of the actual modem config */
        ofdm_config = ofdm_get_config_param();

        ofdm_bitsperframe = ofdm_get_bits_per_frame();
        ofdm_nuwbits = (ofdm_config->ns - 1) * ofdm_config->bps - ofdm_config->txtbits;
        ofdm_ntxtbits = ofdm_config->txtbits;

        f->ldpc = (struct LDPC*)MALLOC(sizeof(struct LDPC));

        if (f->ldpc == NULL) {
            return NULL;
        }

        set_up_hra_112_112(f->ldpc, ofdm_config);
#ifdef __EMBEDDED__
	f->ldpc->max_iter = 10;
#endif	
        int coded_syms_per_frame = f->ldpc->coded_syms_per_frame;
        
        if (adv == NULL) {
            f->interleave_frames = 1;
        } else {
            assert((adv->interleave_frames >= 0) && (adv->interleave_frames <= 16));
            f->interleave_frames = adv->interleave_frames;
        }

        f->modem_frame_count_tx = f->modem_frame_count_rx = 0;
        
        f->codeword_symbols = (COMP*)MALLOC(sizeof(COMP)*f->interleave_frames*coded_syms_per_frame);

        if (f->codeword_symbols == NULL) {return NULL;}

        f->codeword_amps = (float*)MALLOC(sizeof(float)*f->interleave_frames*coded_syms_per_frame);

        if (f->codeword_amps == NULL) {FREE(f->codeword_symbols); return NULL;}

        for (int i=0; i<f->interleave_frames*coded_syms_per_frame; i++) {
            f->codeword_symbols[i].real = 0.0;
            f->codeword_symbols[i].imag = 0.0;
            f->codeword_amps[i] = 0.0;
        }

        f->nin = ofdm_get_samples_per_frame();
        f->n_nat_modem_samples = ofdm_get_samples_per_frame();
        f->n_nom_modem_samples = ofdm_get_samples_per_frame();
        f->n_max_modem_samples = ofdm_get_max_samples_per_frame();
        f->modem_sample_rate = ofdm_config->fs;
        f->clip = 0;

        nbit = f->sz_error_pattern = ofdm_bitsperframe;

        f->tx_bits = NULL; /* not used for 700D */

	if (f->interleave_frames > 1) {
            /* only allocate this array for interleaver sizes > 1 to save memory on SM1000 port */
            f->mod_out = (COMP*)MALLOC(sizeof(COMP)*f->interleave_frames*f->n_nat_modem_samples);
            if (f->mod_out == NULL) {
                if (f->codeword_symbols != NULL) { FREE (f->codeword_symbols); }
                return NULL;
            }

            for (int i=0; i<f->interleave_frames*f->n_nat_modem_samples; i++) {
                f->mod_out[i].real = 0.0;
                f->mod_out[i].imag = 0.0;
            }
	}

#ifndef __EMBEDDED__
        /* tx BPF off on embedded platforms, as it consumes significant CPU */
        ofdm_set_tx_bpf(f->ofdm, 1);
#endif

    }
        
#ifdef __LPCNET__
    if (FDV_MODE_ACTIVE( FREEDV_MODE_2020, mode) ) {
        f->speech_sample_rate = FS_VOICE_16K;
        f->snr_squelch_thresh = 4.0;
        f->squelch_en = 0;
        codec2_mode = CODEC_MODE_LPCNET_1733;
        
        if ((ofdm_config = (struct OFDM_CONFIG *) CALLOC(1, sizeof (struct OFDM_CONFIG))) == NULL) {
            // TODO: do we need this code? Can these actually be allocated before now?
            if (f->tx_bits != NULL) {
              FREE(f->tx_bits);
              f->tx_bits = NULL;
            }
            if (f->rx_bits != NULL) {
              FREE(f->rx_bits);
              f->rx_bits = NULL;
            }
            return NULL;
        }
        
        if (adv == NULL) {
            f->interleave_frames = 1;
        } else {
            assert((adv->interleave_frames >= 0) && (adv->interleave_frames <= 16));
            f->interleave_frames = adv->interleave_frames;
        }

        f->ofdm = ofdm_create(ofdm_config);
        FREE(ofdm_config);

        /* Get a copy of the 700D modem config as template then modify for 2020 */
        ofdm_config = ofdm_get_config_param();
        ofdm_destroy(f->ofdm);
        ofdm_config->nc = 31; int data_bits_per_frame = 312;
        ofdm_config->ts = 0.0205;
        f->ofdm = ofdm_create(ofdm_config);
            
        f->ldpc = (struct LDPC*)MALLOC(sizeof(struct LDPC));
        if (f->ldpc == NULL) {
            return NULL;
        }
        
        set_up_hra_504_396(f->ldpc, ofdm_config);
        set_data_bits_per_frame(f->ldpc, data_bits_per_frame, ofdm_config->bps);
        int coded_syms_per_frame = f->ldpc->coded_syms_per_frame;
        
        // TODO: not sure if these are set up correctly, in Octave we
        // hard code num UW and txt bits, have added asserts to check
        ofdm_bitsperframe = ofdm_get_bits_per_frame();
        ofdm_nuwbits = (ofdm_config->ns - 1) * ofdm_config->bps - ofdm_config->txtbits;
        ofdm_ntxtbits = ofdm_config->txtbits;
        assert(ofdm_nuwbits == 10);
        assert(ofdm_ntxtbits == 4);
        if (f->verbose) {
            fprintf(stderr, "ldpc_data_bits_per_frame = %d\n", f->ldpc->ldpc_data_bits_per_frame);
            fprintf(stderr, "ldpc_coded_bits_per_frame  = %d\n", f->ldpc->ldpc_coded_bits_per_frame);
            fprintf(stderr, "data_bits_per_frame = %d\n", data_bits_per_frame);
            fprintf(stderr, "coded_bits_per_frame  = %d\n", f->ldpc->coded_bits_per_frame);
            fprintf(stderr, "coded_syms_per_frame  = %d\n", f->ldpc->coded_syms_per_frame);
            fprintf(stderr, "ofdm_bits_per_frame  = %d\n", ofdm_bitsperframe);
            fprintf(stderr, "interleave_frames: %d\n", f->interleave_frames);
        }
        
        if (adv == NULL) {
            f->interleave_frames = 1;
        } else {
            assert((adv->interleave_frames >= 0) && (adv->interleave_frames <= 16));
            f->interleave_frames = adv->interleave_frames;
        }
        
        f->modem_frame_count_tx = f->modem_frame_count_rx = 0;
        
        f->codeword_symbols = (COMP*)MALLOC(sizeof(COMP)*f->interleave_frames*coded_syms_per_frame);
        if (f->codeword_symbols == NULL) {return NULL;}

        f->codeword_amps = (float*)MALLOC(sizeof(float)*f->interleave_frames*coded_syms_per_frame);
        if (f->codeword_amps == NULL) {FREE(f->codeword_symbols); return NULL;}

        for (int i=0; i<f->interleave_frames*coded_syms_per_frame; i++) {
            f->codeword_symbols[i].real = 0.0;
            f->codeword_symbols[i].imag = 0.0;
            f->codeword_amps[i] = 0.0;
        }

        f->nin = ofdm_get_samples_per_frame();
        f->n_nat_modem_samples = ofdm_get_samples_per_frame();
        f->n_nom_modem_samples = ofdm_get_samples_per_frame();
        f->n_max_modem_samples = ofdm_get_max_samples_per_frame();
        f->modem_sample_rate = ofdm_config->fs;
        f->clip = 0;

        nbit = f->sz_error_pattern = ofdm_bitsperframe;

        f->tx_bits = NULL; /* not used for 2020 */

	if (f->interleave_frames > 1) {
            /* only allocate this array for interleaver sizes > 1 to save memory on SM1000 port */
            f->mod_out = (COMP*)MALLOC(sizeof(COMP)*f->interleave_frames*f->n_nat_modem_samples);
            if (f->mod_out == NULL) {
                if (f->codeword_symbols != NULL) { FREE (f->codeword_symbols); }
                return NULL;
            }

            for (int i=0; i<f->interleave_frames*f->n_nat_modem_samples; i++) {
                f->mod_out[i].real = 0.0;
                f->mod_out[i].imag = 0.0;
            }
	}
        
        /* TODO: tx BPF off by default, as we need new filter coeffs for FreeDV 2020 waveform */
        ofdm_set_tx_bpf(f->ofdm, 0);
        
        codec2_mode = CODEC_MODE_LPCNET_1733; // not really codec 2 but got to call it something        
    }

    /* Malloc something to appease freedv_destroy */
    f->codec_bits = NULL;
#endif
    
    if (FDV_MODE_ACTIVE( FREEDV_MODE_2400A, mode) || FDV_MODE_ACTIVE( FREEDV_MODE_2400B, mode)) {
      
        /* Set up the C2 mode */
        codec2_mode = CODEC2_MODE_1300;
        /* Set the number of protocol bits */
        f->n_protocol_bits = 20;
        f->sz_error_pattern = 0;
        f->ext_vco = 0;
    }
    
    if (FDV_MODE_ACTIVE( FREEDV_MODE_2400A, mode)) {
        /* Create the framer|deframer */
        f->deframer = fvhff_create_deframer(FREEDV_VHF_FRAME_A,0);
        if(f->deframer == NULL)
            return NULL;
  
        f->fsk = fsk_create_hbr(48000,1200,10,4,1200,1200);
        
        /* Note: fsk expects tx/rx bits as an array of uint8_ts, not ints */
        f->tx_bits = (int*)MALLOC(f->fsk->Nbits*sizeof(uint8_t));
        
        if(f->fsk == NULL){
            fvhff_destroy_deframer(f->deframer);
            return NULL;
        }
        
        f->n_nom_modem_samples = f->fsk->N;
        f->n_max_modem_samples = f->fsk->N + (f->fsk->Ts);
        f->n_nat_modem_samples = f->fsk->N;
        f->nin = fsk_nin(f->fsk);
        f->modem_sample_rate = 48000;
        f->modem_symbol_rate = 1200;
        /* Malloc something to appease freedv_init and freedv_destroy */
        f->codec_bits = MALLOC(1);
    }
    
    if (FDV_MODE_ACTIVE( FREEDV_MODE_2400B, mode) ) {
        /* Create the framer|deframer */
        f->deframer = fvhff_create_deframer(FREEDV_VHF_FRAME_A,1);
        if(f->deframer == NULL) {
            if (f->codec_bits != NULL) { FREE(f->codec_bits); }
            return NULL;
        }
        
        f->fmfsk = fmfsk_create(48000,2400);
         
        if(f->fmfsk == NULL){
            fvhff_destroy_deframer(f->deframer);
            return NULL;
        }
        /* Note: fsk expects tx/rx bits as an array of uint8_ts, not ints */
        f->tx_bits = (int*)MALLOC(f->fmfsk->nbit*sizeof(uint8_t));
        
        f->n_nom_modem_samples = f->fmfsk->N;
        f->n_max_modem_samples = f->fmfsk->N + (f->fmfsk->Ts);
        f->n_nat_modem_samples = f->fmfsk->N;
        f->nin = fmfsk_nin(f->fmfsk);
        f->modem_sample_rate = 48000;
        /* Malloc something to appease freedv_init and freedv_destroy */
        f->codec_bits = MALLOC(1);
    }
    
    if (FDV_MODE_ACTIVE( FREEDV_MODE_800XA, mode)) {
        /* Create the framer|deframer */
        f->deframer = fvhff_create_deframer(FREEDV_HF_FRAME_B,0);
        if(f->deframer == NULL)
            return NULL;
  
        f->fsk = fsk_create_hbr(8000,400,10,4,800,400);
        fsk_set_nsym(f->fsk,32);
        
        /* Note: fsk expects tx/rx bits as an array of uint8_ts, not ints */
        f->tx_bits = (int*)MALLOC(f->fsk->Nbits*sizeof(uint8_t));
        
        if(f->fsk == NULL){
            fvhff_destroy_deframer(f->deframer);
            return NULL;
        }
        
        f->n_nom_modem_samples = f->fsk->N;
        f->n_max_modem_samples = f->fsk->N + (f->fsk->Ts);
        f->n_nat_modem_samples = f->fsk->N;
        f->nin = fsk_nin(f->fsk);
        f->modem_sample_rate = 8000;
        f->modem_symbol_rate = 400;

        /* Malloc something to appease freedv_init and freedv_destroy */
        f->codec_bits = MALLOC(1);
        
        f->n_protocol_bits = 0;
        codec2_mode = CODEC2_MODE_700C;
        fsk_stats_normalise_eye(f->fsk, 0);
        f->sz_error_pattern = 0;
    }

    /* Init test frame states */
    
    f->test_frames_diversity = 1;
    f->test_frame_sync_state = 0;
    f->test_frame_sync_state_upper = 0;
    f->total_bits = 0;
    f->total_bit_errors = 0;
    f->total_bits_coded = 0;
    f->total_bit_errors_coded = 0;

    /* -----------------------------------------------------------------------------------------------*\
    |                         Init Codec for this FreeDV mode                                          |
    \*------------------------------------------------------------------------------------------------*/

    if (codec2_mode == CODEC_MODE_LPCNET_1733) {
#ifdef __LPCNET__
        f->lpcnet = lpcnet_freedv_create(1);
#endif
        f->codec2 = NULL;
    }
    else {
        f->codec2 = codec2_create(codec2_mode);
        if (f->codec2 == NULL)
            return NULL;
    }
    
    /* work out how many codec 2 frames per mode frame, and number of
       bytes of storage for packed and unpacked bits.  TODO: do we really
       need to work in packed bits at all?  It's messy, chars would probably
       be OK.... */
    
    if (FDV_MODE_ACTIVE( FREEDV_MODE_1600, mode) || FDV_MODE_ACTIVE( FREEDV_MODE_2400A, mode) || FDV_MODE_ACTIVE( FREEDV_MODE_2400B, mode)) {
        f->n_speech_samples = codec2_samples_per_frame(f->codec2);
        f->n_codec_bits = codec2_bits_per_frame(f->codec2);
        nbit = f->n_codec_bits;
        nbyte = (nbit + 7) / 8;
    } else if (FDV_MODE_ACTIVE( FREEDV_MODE_800XA, mode) ) {
        f->n_speech_samples = 2*codec2_samples_per_frame(f->codec2);
        f->n_codec_bits = codec2_bits_per_frame(f->codec2);
        nbit = f->n_codec_bits;
        nbyte = (nbit + 7) / 8;
        nbyte = nbyte*2;
        nbit = 8*nbyte;
        f->n_codec_bits = nbit;
    } else if (FDV_MODE_ACTIVE( FREEDV_MODE_700, mode) || FDV_MODE_ACTIVE( FREEDV_MODE_700B, mode) || FDV_MODE_ACTIVE( FREEDV_MODE_700C, mode)) {
        f->n_speech_samples = 2*codec2_samples_per_frame(f->codec2);
        f->n_codec_bits = 2*codec2_bits_per_frame(f->codec2);
        nbit = f->n_codec_bits;
        nbyte = 2*((codec2_bits_per_frame(f->codec2) + 7) / 8);
    } else if (FDV_MODE_ACTIVE(FREEDV_MODE_700D, mode)) /* mode == FREEDV_MODE_700D */ {

        /* should be exactly an integer number of Codec 2 frames in a OFDM modem frame */

        assert((f->ldpc->data_bits_per_frame % codec2_bits_per_frame(f->codec2)) == 0);

        int Ncodec2frames = f->ldpc->data_bits_per_frame/codec2_bits_per_frame(f->codec2);

        f->n_speech_samples = Ncodec2frames*codec2_samples_per_frame(f->codec2);
        f->n_codec_bits = f->interleave_frames*Ncodec2frames*codec2_bits_per_frame(f->codec2);
        nbit = codec2_bits_per_frame(f->codec2);
        nbyte = (nbit + 7) / 8;
        nbyte = nbyte*Ncodec2frames*f->interleave_frames;
        f->nbyte_packed_codec_bits = nbyte;

        //fprintf(stderr, "Ncodec2frames: %d n_speech_samples: %d n_codec_bits: %d nbit: %d  nbyte: %d\n",
        //        Ncodec2frames, f->n_speech_samples, f->n_codec_bits, nbit, nbyte);

        f->packed_codec_bits_tx = (unsigned char*)MALLOC(nbyte*sizeof(char));
        if (f->packed_codec_bits_tx == NULL) return(NULL);
        f->codec_bits = NULL;
    } else if (FDV_MODE_ACTIVE(FREEDV_MODE_2020, mode)) {
#ifdef __LPCNET__
        /* should be exactly an integer number of Codec frames in a OFDM modem frame */

        assert((f->ldpc->data_bits_per_frame % lpcnet_bits_per_frame(f->lpcnet)) == 0);

        int Ncodecframes = f->ldpc->data_bits_per_frame/lpcnet_bits_per_frame(f->lpcnet);

        f->n_speech_samples = Ncodecframes*lpcnet_samples_per_frame(f->lpcnet);
        f->n_codec_bits = f->interleave_frames*Ncodecframes*lpcnet_bits_per_frame(f->lpcnet);
        nbit = f->n_codec_bits;

        // we actually have unpacked data but uses this as it's convenient
        nbyte = nbit;
#endif
    }
    
    f->packed_codec_bits = (unsigned char*)MALLOC(nbyte*sizeof(char));
    if (f->packed_codec_bits == NULL) return(NULL);

    if (FDV_MODE_ACTIVE( FREEDV_MODE_1600, mode))
        f->codec_bits = (int*)MALLOC(nbit*sizeof(int));

    if (FDV_MODE_ACTIVE( FREEDV_MODE_700, mode) || FDV_MODE_ACTIVE( FREEDV_MODE_700B, mode) || FDV_MODE_ACTIVE( FREEDV_MODE_700C, mode))
        f->codec_bits = (int*)MALLOC(COHPSK_BITS_PER_FRAME*sizeof(int));
        
    /* Note: VHF Framer/deframer goes directly from packed codec/vc/proto bits to filled frame */

    if (f->packed_codec_bits == NULL) {
        if (f->codec_bits != NULL) {FREE(f->codec_bits); f->codec_bits = NULL; }
        return NULL;
    }

    /* Sample rate conversion for modes using COHPSK */
    
    if (FDV_MODE_ACTIVE( FREEDV_MODE_700, mode) || FDV_MODE_ACTIVE( FREEDV_MODE_700B, mode) || FDV_MODE_ACTIVE( FREEDV_MODE_700C, mode) ) { 
        f->ptFilter7500to8000 = (struct quisk_cfFilter *)MALLOC(sizeof(struct quisk_cfFilter));
        f->ptFilter8000to7500 = (struct quisk_cfFilter *)MALLOC(sizeof(struct quisk_cfFilter));
        quisk_filt_cfInit(f->ptFilter8000to7500, quiskFilt120t480, sizeof(quiskFilt120t480)/sizeof(float));
        quisk_filt_cfInit(f->ptFilter7500to8000, quiskFilt120t480, sizeof(quiskFilt120t480)/sizeof(float));
    } else {
        f->ptFilter7500to8000 = NULL;
        f->ptFilter8000to7500 = NULL;
    }

    /* Varicode low bit rate text states */
    
    varicode_decode_init(&f->varicode_dec_states, 1);
    f->nvaricode_bits = 0;
    f->varicode_bit_index = 0;
    f->freedv_get_next_tx_char = NULL;
    f->freedv_put_next_rx_char = NULL;
    f->freedv_put_next_proto = NULL;
    f->freedv_get_next_proto = NULL;
    f->total_bit_errors = 0;

    return f;
}

/*---------------------------------------------------------------------------*\

  FUNCTION....: freedv_close
  AUTHOR......: David Rowe
  DATE CREATED: 3 August 2014

  Frees up memory.

\*---------------------------------------------------------------------------*/

void freedv_close(struct freedv *freedv) {
    assert(freedv != NULL);

    FREE(freedv->packed_codec_bits);
    FREE(freedv->codec_bits);
    FREE(freedv->tx_bits);
    FREE(freedv->rx_bits);
    if (FDV_MODE_ACTIVE(FREEDV_MODE_1600, freedv->mode)) {
        FREE(freedv->fdmdv_bits);
        fdmdv_destroy(freedv->fdmdv);
    }
    if (FDV_MODE_ACTIVE( FREEDV_MODE_700, freedv->mode) || FDV_MODE_ACTIVE( FREEDV_MODE_700B, freedv->mode) || FDV_MODE_ACTIVE( FREEDV_MODE_700C, freedv->mode))
        cohpsk_destroy(freedv->cohpsk);
    if (FDV_MODE_ACTIVE( FREEDV_MODE_700D, freedv->mode)) {
        FREE(freedv->packed_codec_bits_tx);
        if (freedv->interleave_frames > 1)
            FREE(freedv->mod_out);
        FREE(freedv->codeword_symbols);
        FREE(freedv->codeword_amps);
        FREE(freedv->ldpc);
        ofdm_destroy(freedv->ofdm);
    }
    if (FDV_MODE_ACTIVE( FREEDV_MODE_2020, freedv->mode)) {
        if (freedv->interleave_frames > 1)
            FREE(freedv->mod_out);
        FREE(freedv->codeword_symbols);
        FREE(freedv->codeword_amps);
        FREE(freedv->ldpc);
        ofdm_destroy(freedv->ofdm);
#ifdef __LPCNET__
        lpcnet_freedv_destroy(freedv->lpcnet);
#endif        
    }
    if (FDV_MODE_ACTIVE( FREEDV_MODE_2400A, freedv->mode) || FDV_MODE_ACTIVE( FREEDV_MODE_800XA, freedv->mode)){
        fsk_destroy(freedv->fsk);
        fvhff_destroy_deframer(freedv->deframer);
    }
    
    if (FDV_MODE_ACTIVE( FREEDV_MODE_2400B, freedv->mode)){
        fmfsk_destroy(freedv->fmfsk);
		fvhff_destroy_deframer(freedv->deframer);
    }
    
    if (freedv->codec2)
        codec2_destroy(freedv->codec2);
    if (freedv->ptFilter8000to7500) {
        quisk_filt_destroy(freedv->ptFilter8000to7500);
        FREE(freedv->ptFilter8000to7500);
        freedv->ptFilter8000to7500 = NULL;
    }
    if (freedv->ptFilter7500to8000) {
        quisk_filt_destroy(freedv->ptFilter7500to8000);
        FREE(freedv->ptFilter7500to8000);
        freedv->ptFilter7500to8000 = NULL;
    }
    FREE(freedv);
}

/*---------------------------------------------------------------------------*\

  FUNCTION....: freedv_tx
  AUTHOR......: David Rowe
  DATE CREATED: 3 August 2014

  Takes a frame of input speech samples, encodes and modulates them to
  produce a frame of modem samples that can be sent to the
  transmitter.  See freedv_tx.c for an example.

  speech_in[] is sampled at 8000 Hz, and the user must supply a block
  of exactly freedv_get_n_speech_samples(). The speech_in[] level
  should be such that the peak speech level is between +/- 16384 and
  +/- 32767.

  The FDM modem signal mod_out[] is sampled at 8000 Hz and is
  freedv_get_n_nom_modem_samples() long.  mod_out[] will be scaled
  such that the peak level is just less than +/-32767.

  The complex-valued output can directly drive an I/Q modulator to
  produce a single sideband signal.  To generate the other sideband,
  take the complex conjugate of mod_out[].

  The FreeDV 1600 modem has a high crest factor (around 12dB), however
  the energy and duration of the peaks is small.  FreeDV 1600 is
  usually operated at a "backoff" of 8dB.  Adjust the power amplifier
  drive so that the average power is 8dB less than the peak power of
  the PA.  For example, on a radio rated at 100W PEP for SSB, the
  average FreeDV power is typically 20W.

  The FreeDV 700 modem has a crest factor of about 8dB (with
  f->clip=1, the default), so if your PA can handle it, it can be
  driven harder than FreeDV 1600.  Caution - some PAs cannot handle a
  high continuous power.  A conservative level is 20W average for a
  100W PEP rated PA.

\*---------------------------------------------------------------------------*/

/* real-valued short sample output, useful for going straight to DAC */

/* TX routines for 2400 FSK modes, after codec2 encoding */
static void freedv_tx_fsk_voice(struct freedv *f, short mod_out[]) {
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
            fvhff_frame_bits(FREEDV_VHF_FRAME_A,(uint8_t*)(f->tx_bits),(uint8_t*)(f->packed_codec_bits),proto_bits,vc_bits);
        }else if(f->freedv_get_next_tx_char != NULL){
            fvhff_frame_bits(FREEDV_VHF_FRAME_A,(uint8_t*)(f->tx_bits),(uint8_t*)(f->packed_codec_bits),NULL,vc_bits);
        }else {
            fvhff_frame_bits(FREEDV_VHF_FRAME_A,(uint8_t*)(f->tx_bits),(uint8_t*)(f->packed_codec_bits),NULL,NULL);
        }
    /* Frame for 800XA */
    }else if(FDV_MODE_ACTIVE( FREEDV_MODE_800XA, f->mode)){
        fvhff_frame_bits(FREEDV_HF_FRAME_B,(uint8_t*)(f->tx_bits),(uint8_t*)(f->packed_codec_bits),NULL,NULL);
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
static void freedv_comptx_fsk_voice(struct freedv *f, COMP mod_out[]) {
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
            fvhff_frame_bits(FREEDV_VHF_FRAME_A,(uint8_t*)(f->tx_bits),(uint8_t*)(f->packed_codec_bits),proto_bits,vc_bits);
        }else if(f->freedv_get_next_tx_char != NULL){
            fvhff_frame_bits(FREEDV_VHF_FRAME_A,(uint8_t*)(f->tx_bits),(uint8_t*)(f->packed_codec_bits),NULL,vc_bits);
        }else {
            fvhff_frame_bits(FREEDV_VHF_FRAME_A,(uint8_t*)(f->tx_bits),(uint8_t*)(f->packed_codec_bits),NULL,NULL);
        }
    /* Frame for 800XA */
    }else if(FDV_MODE_ACTIVE( FREEDV_MODE_800XA, f->mode)){
        fvhff_frame_bits(FREEDV_HF_FRAME_B,(uint8_t*)(f->tx_bits),(uint8_t*)(f->packed_codec_bits),NULL,NULL);
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
static void freedv_tx_fsk_data(struct freedv *f, short mod_out[]) {
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

void freedv_tx(struct freedv *f, short mod_out[], short speech_in[]) {
    assert(f != NULL);
    COMP tx_fdm[f->n_nom_modem_samples];
    int  i;
    assert((FDV_MODE_ACTIVE( FREEDV_MODE_1600, f->mode))  || (FDV_MODE_ACTIVE( FREEDV_MODE_700, f->mode))   || 
           (FDV_MODE_ACTIVE( FREEDV_MODE_700B, f->mode))  || (FDV_MODE_ACTIVE( FREEDV_MODE_700C, f->mode))  ||
           (FDV_MODE_ACTIVE( FREEDV_MODE_700D, f->mode))  || 
           (FDV_MODE_ACTIVE( FREEDV_MODE_2400A, f->mode)) || (FDV_MODE_ACTIVE( FREEDV_MODE_2400B, f->mode)) || 
           (FDV_MODE_ACTIVE( FREEDV_MODE_800XA, f->mode)) || (FDV_MODE_ACTIVE( FREEDV_MODE_2020, f->mode)));
    
    /* FSK and MEFSK/FMFSK modems work only on real samples. It's simpler to just 
     * stick them in the real sample tx/rx functions than to add a comp->real converter
     * to comptx */
     
    if ((FDV_MODE_ACTIVE( FREEDV_MODE_2400A, f->mode)) || (FDV_MODE_ACTIVE( FREEDV_MODE_2400B, f->mode)) || (FDV_MODE_ACTIVE( FREEDV_MODE_800XA, f->mode))){
        /* 800XA has two codec frames per modem frame */
        if(FDV_MODE_ACTIVE( FREEDV_MODE_800XA, f->mode)){
            codec2_encode(f->codec2, &f->packed_codec_bits[0], &speech_in[  0]);
            codec2_encode(f->codec2, &f->packed_codec_bits[4], &speech_in[320]);
        }else{
            codec2_encode(f->codec2, f->packed_codec_bits, speech_in);
        }
        freedv_tx_fsk_voice(f, mod_out);
    } else{
        freedv_comptx(f, tx_fdm, speech_in);
        for(i=0; i<f->n_nom_modem_samples; i++)
            mod_out[i] = tx_fdm[i].real;
    }
}

/* complex valued output, useful for suitable for single sided freq shifting */

static void freedv_comptx_fdmdv_1600(struct freedv *f, COMP mod_out[]) {
    int    bit, byte, i, j;
    int    bits_per_codec_frame, bits_per_modem_frame;
    int    data, codeword1, data_flag_index;
    COMP   tx_fdm[f->n_nat_modem_samples];

    bits_per_codec_frame = codec2_bits_per_frame(f->codec2);
    bits_per_modem_frame = fdmdv_bits_per_frame(f->fdmdv);

    /* unpack bits, MSB first */

    bit = 7; byte = 0;
    for(i=0; i<bits_per_codec_frame; i++) {
        f->codec_bits[i] = (f->packed_codec_bits[byte] >> bit) & 0x1;
        bit--;
        if (bit < 0) {
            bit = 7;
            byte++;
        }
    }

    // spare bit in frame that codec defines.  Use this 1
    // bit/frame to send txt messages

    data_flag_index = codec2_get_spare_bit_index(f->codec2);

    if (f->nvaricode_bits) {
        f->codec_bits[data_flag_index] = f->tx_varicode_bits[f->varicode_bit_index++];
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

    /* Protect first 12 out of first 16 excitation bits with (23,12) Golay Code:

       0,1,2,3: v[0]..v[1]
       4,5,6,7: MSB of pitch
       11,12,13,14: MSB of energy

    */

    data = 0;
    for(i=0; i<8; i++) {
        data <<= 1;
        data |= f->codec_bits[i];
    }
    for(i=11; i<15; i++) {
        data <<= 1;
        data |= f->codec_bits[i];
    }
    codeword1 = golay23_encode(data);

    /* now pack output frame with parity bits at end to make them
       as far apart as possible from the data they protect.  Parity
       bits are LSB of the Golay codeword */

    for(i=0; i<bits_per_codec_frame; i++)
        f->tx_bits[i] = f->codec_bits[i];
    for(j=0; i<bits_per_codec_frame+11; i++,j++) {
        f->tx_bits[i] = (codeword1 >> (10-j)) & 0x1;
    }
    f->tx_bits[i] = 0; /* spare bit */

    /* optionally overwrite with test frames */

    if (f->test_frames) {
        fdmdv_get_test_bits(f->fdmdv, f->tx_bits);
        fdmdv_get_test_bits(f->fdmdv, &f->tx_bits[bits_per_modem_frame]);
        //fprintf(stderr, "test frames on tx\n");
    }

    /* modulate even and odd frames */

    fdmdv_mod(f->fdmdv, tx_fdm, f->tx_bits, &f->tx_sync_bit);
    assert(f->tx_sync_bit == 1);

    fdmdv_mod(f->fdmdv, &tx_fdm[FDMDV_NOM_SAMPLES_PER_FRAME], &f->tx_bits[bits_per_modem_frame], &f->tx_sync_bit);
    assert(f->tx_sync_bit == 0);

    assert(2*FDMDV_NOM_SAMPLES_PER_FRAME == f->n_nom_modem_samples);

    for(i=0; i<f->n_nom_modem_samples; i++)
        mod_out[i] = fcmult(FDMDV_SCALE, tx_fdm[i]);
}

static void freedv_comptx_700(struct freedv *f, COMP mod_out[]) {
    int    bit, byte, i, j, k;
    int    bits_per_codec_frame, bits_per_modem_frame;
    int    data_flag_index, nspare;
    COMP   tx_fdm[f->n_nat_modem_samples];

    bits_per_codec_frame = codec2_bits_per_frame(f->codec2);
    bits_per_modem_frame = COHPSK_BITS_PER_FRAME;

    byte = 0;
    for (j=0; j<bits_per_modem_frame; j+=bits_per_codec_frame) {

        /* unpack bits, MSB first */

        bit = 7;
        for(i=0; i<bits_per_codec_frame; i++) {
            f->codec_bits[j+i] = (f->packed_codec_bits[byte] >> bit) & 0x1;
            bit--;
            if (bit < 0) {
                bit = 7;
                byte++;
            }
        }
	if (bit != 7)
	    byte++;

        // spare bits in frame that codec defines.  Use these spare
        // bits/frame to send txt messages

        switch(f->mode) {
        case FREEDV_MODE_700:
            nspare = 2;
            break;
        case FREEDV_MODE_700B:
            nspare = 1; // Just one spare bit for FREEDV_MODE_700B
            break;
        case FREEDV_MODE_700C:
            nspare = 0; // and no spare bits for 700C atm
            break;
        default:
            assert(0);
        }

        data_flag_index = codec2_get_spare_bit_index(f->codec2);

        for(k=0; k<nspare; k++) {
            if (f->nvaricode_bits) {
                f->codec_bits[j+data_flag_index+k] = f->tx_varicode_bits[f->varicode_bit_index++];
                //fprintf(stderr, "%d %d\n", j+data_flag_index+k, f->codec_bits[j+data_flag_index+k]);
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
    }
    /* optionally ovwerwrite the codec bits with test frames */

    if (f->test_frames) {
        cohpsk_get_test_bits(f->cohpsk, f->codec_bits);
    }

    /* cohpsk modulator */

    cohpsk_mod(f->cohpsk, tx_fdm, f->codec_bits, COHPSK_BITS_PER_FRAME);
    if (f->clip) 
        cohpsk_clip(tx_fdm, COHPSK_CLIP, COHPSK_NOM_SAMPLES_PER_FRAME);
    for(i=0; i<f->n_nat_modem_samples; i++)
        mod_out[i] = fcmult(FDMDV_SCALE*NORM_PWR_COHPSK, tx_fdm[i]);
    i = quisk_cfInterpDecim((complex float *)mod_out, f->n_nat_modem_samples, f->ptFilter7500to8000, 16, 15);
    //assert(i == f->n_nom_modem_samples);
    // Caution: assert fails if f->n_nat_modem_samples * 16.0 / 15.0 is not an integer

}

/*
  Ok so when interleaved, we take the interleaver length of input samples,
  and output that many modem samples, e.g. for interleaver of length 4:

  record input speech 1234
  freedv tx              |
  play modem sig         1234
  record modem sig       1234
  freedv_rx                   |
  play output speech          1234        
  time axis --------->123456789012---->

  So a sample of input speech at time 1 is ouput at time 9.  We assume
  the freedv_tx and freedv_rx and propogation time over channel (from
  when a modem signal is played at the HF tx to when it is recorded at
  the HF rx) is zero.

  The freedv tx interface ouputs n_nom_modem_samples, which a single
  OFDM modem frame, 112 payload bits or 4 speech codec frames.  So
  this function must always have 1280 speech samples as input, and
  1280 modem samples as output, regradless of interleaver_frames.  For
  interleaver_frames > 1, we need to buffer samples.
*/

static void freedv_comptx_700d(struct freedv *f, COMP mod_out[]) {
    int    bit, byte, i, j, k;
    int    nspare;
 
    int data_bits_per_frame = f->ldpc->data_bits_per_frame;
    int bits_per_interleaved_frame = f->interleave_frames*data_bits_per_frame;
    uint8_t tx_bits[bits_per_interleaved_frame];
    int bits_per_codec_frame = codec2_bits_per_frame(f->codec2);

    byte = 0;
    for (j=0; j<bits_per_interleaved_frame; j+=bits_per_codec_frame) {

        /* unpack bits, MSB first */

        bit = 7;
        for(i=0; i<bits_per_codec_frame; i++) {
            tx_bits[j+i] = (f->packed_codec_bits_tx[byte] >> bit) & 0x1;
            bit--;
            if (bit < 0) {
                bit = 7;
                byte++;
            }
        }
	if (bit != 7)
	    byte++;
    }

    assert(byte <= f->nbyte_packed_codec_bits);
    
    // Generate Varicode txt bits. Txt bits in OFDM frame come just
    // after Unique Word (UW).  Txt bits aren't protected by FEC, and need to be
    // added to each frame after interleaver as done it's thing

    nspare = ofdm_ntxtbits*f->interleave_frames;
    uint8_t txt_bits[nspare];

    for(k=0; k<nspare; k++) {
        if (f->nvaricode_bits == 0) {
            /* get new char and encode */
            char s[2];
            if (f->freedv_get_next_tx_char != NULL) {
                s[0] = (*f->freedv_get_next_tx_char)(f->callback_state);
                f->nvaricode_bits = varicode_encode(f->tx_varicode_bits, s, VARICODE_MAX_BITS, 1, 1);
                f->varicode_bit_index = 0;
            }
        }
        if (f->nvaricode_bits) {
            txt_bits[k] = f->tx_varicode_bits[f->varicode_bit_index++];
            f->nvaricode_bits--;
        }
    }

    /* optionally replace codec payload bits with test frames known to rx */

    if (f->test_frames) {
        uint8_t payload_data_bits[data_bits_per_frame];
        ofdm_generate_payload_data_bits(payload_data_bits, data_bits_per_frame);

        for (j=0; j<f->interleave_frames; j++) {
            for(i=0; i<data_bits_per_frame; i++) {
                tx_bits[j*data_bits_per_frame + i] = payload_data_bits[i];
            }
        }
    }

    /* OK now ready to LDPC encode, interleave, and OFDM modulate */
    
    complex float tx_sams[f->interleave_frames*f->n_nat_modem_samples];
    COMP asam;
    
    ofdm_ldpc_interleave_tx(f->ofdm, f->ldpc, tx_sams, tx_bits, txt_bits, f->interleave_frames, ofdm_config);

    for(i=0; i<f->interleave_frames*f->n_nat_modem_samples; i++) {
        asam.real = crealf(tx_sams[i]);
        asam.imag = cimagf(tx_sams[i]);
        mod_out[i] = fcmult(OFDM_AMP_SCALE*NORM_PWR_OFDM, asam);
    }

    if (f->clip) {
        //fprintf(stderr, "clip ");
        cohpsk_clip(mod_out, OFDM_CLIP, f->interleave_frames*f->n_nat_modem_samples);
    }
}


#ifdef __LPCNET__
static void freedv_comptx_2020(struct freedv *f, COMP mod_out[]) {
    int    i, j, k;
    int    nspare;
 
    int data_bits_per_frame = f->ldpc->data_bits_per_frame;
    int bits_per_interleaved_frame = f->interleave_frames*data_bits_per_frame;
    uint8_t tx_bits[bits_per_interleaved_frame];

    memcpy(tx_bits, f->packed_codec_bits, bits_per_interleaved_frame);
    
    // Generate Varicode txt bits. Txt bits in OFDM frame come just
    // after Unique Word (UW).  Txt bits aren't protected by FEC, and need to be
    // added to each frame after interleaver as done it's thing

    nspare = ofdm_ntxtbits*f->interleave_frames;
    uint8_t txt_bits[nspare];

    for(k=0; k<nspare; k++) {
        if (f->nvaricode_bits == 0) {
            /* get new char and encode */
            char s[2];
            if (f->freedv_get_next_tx_char != NULL) {
                s[0] = (*f->freedv_get_next_tx_char)(f->callback_state);
                f->nvaricode_bits = varicode_encode(f->tx_varicode_bits, s, VARICODE_MAX_BITS, 1, 1);
                f->varicode_bit_index = 0;
            }
        }
        if (f->nvaricode_bits) {
            txt_bits[k] = f->tx_varicode_bits[f->varicode_bit_index++];
            f->nvaricode_bits--;
        }
    }

    /* optionally replace codec payload bits with test frames known to rx */

    if (f->test_frames) {
        uint8_t payload_data_bits[data_bits_per_frame];
        ofdm_generate_payload_data_bits(payload_data_bits, data_bits_per_frame);

        for (j=0; j<f->interleave_frames; j++) {
            for(i=0; i<data_bits_per_frame; i++) {
                tx_bits[j*data_bits_per_frame + i] = payload_data_bits[i];
            }
        }
    }

    /* OK now ready to LDPC encode, interleave, and OFDM modulate */
    
    complex float tx_sams[f->interleave_frames*f->n_nat_modem_samples];
    COMP asam;
    
    ofdm_ldpc_interleave_tx(f->ofdm, f->ldpc, tx_sams, tx_bits, txt_bits, f->interleave_frames, ofdm_config);

    for(i=0; i<f->interleave_frames*f->n_nat_modem_samples; i++) {
        asam.real = crealf(tx_sams[i]);
        asam.imag = cimagf(tx_sams[i]);
        mod_out[i] = fcmult(OFDM_AMP_SCALE*NORM_PWR_OFDM, asam);
    }

    if (f->clip) {
        cohpsk_clip(mod_out, OFDM_CLIP, f->interleave_frames*f->n_nat_modem_samples);
    }
}
#endif


void freedv_comptx(struct freedv *f, COMP mod_out[], short speech_in[]) {
    assert(f != NULL);

    assert((FDV_MODE_ACTIVE( FREEDV_MODE_1600, f->mode)) || (FDV_MODE_ACTIVE( FREEDV_MODE_700, f->mode)) || 
           (FDV_MODE_ACTIVE( FREEDV_MODE_700B, f->mode)) || (FDV_MODE_ACTIVE( FREEDV_MODE_700C, f->mode)) || 
           (FDV_MODE_ACTIVE( FREEDV_MODE_2400A, f->mode)) || (FDV_MODE_ACTIVE( FREEDV_MODE_2400B, f->mode)) ||
           (FDV_MODE_ACTIVE( FREEDV_MODE_700D, f->mode))  || (FDV_MODE_ACTIVE( FREEDV_MODE_2020, f->mode)));

    if (FDV_MODE_ACTIVE( FREEDV_MODE_1600, f->mode)) {
        codec2_encode(f->codec2, f->packed_codec_bits, speech_in);
        freedv_comptx_fdmdv_1600(f, mod_out);
    }

    int bits_per_codec_frame=0; int bytes_per_codec_frame=0;
    if (f->codec2) {
        bits_per_codec_frame = codec2_bits_per_frame(f->codec2);
        bytes_per_codec_frame = (bits_per_codec_frame + 7) / 8;
    }
    int i,j;
    
    /* all these modes need to pack a bunch of codec frames into one modem frame */
    
    if ((FDV_MODE_ACTIVE( FREEDV_MODE_700, f->mode)) || (FDV_MODE_ACTIVE( FREEDV_MODE_700B, f->mode)) || (FDV_MODE_ACTIVE( FREEDV_MODE_700C, f->mode))) {
	int codec_frames = f->n_codec_bits / bits_per_codec_frame;

        for (j=0; j<codec_frames; j++) {
            codec2_encode(f->codec2, f->packed_codec_bits + j * bytes_per_codec_frame, speech_in);
            speech_in += codec2_samples_per_frame(f->codec2);
        }
        freedv_comptx_700(f, mod_out);
    }

    /* special treatment due to interleaver */
    
    if (FDV_MODE_ACTIVE( FREEDV_MODE_700D, f->mode)) {
        int data_bits_per_frame = f->ldpc->data_bits_per_frame;
	int codec_frames = data_bits_per_frame / bits_per_codec_frame;

        //fprintf(stderr, "modem_frame_count_tx: %d dec_frames: %d bytes offset: %d\n",
        //        f->modem_frame_count_tx, codec_frames, (f->modem_frame_count_tx*codec_frames)*bytes_per_codec_frame);
       
        /* buffer up bits until we get enough encoded bits for interleaver */
        
        for (j=0; j<codec_frames; j++) {
            codec2_encode(f->codec2, f->packed_codec_bits_tx + (f->modem_frame_count_tx*codec_frames+j)*bytes_per_codec_frame, speech_in);
            speech_in += codec2_samples_per_frame(f->codec2);
        }


	/* Only use extra local buffer if needed for interleave > 1 */
	if (f->interleave_frames == 1) {
            freedv_comptx_700d(f, mod_out);
	} else {
            /* call modulate function when we have enough frames to run interleaver */
            assert((f->modem_frame_count_tx >= 0) && 
	    		(f->modem_frame_count_tx < f->interleave_frames));
            f->modem_frame_count_tx++;
            if (f->modem_frame_count_tx == f->interleave_frames) {
                freedv_comptx_700d(f, f->mod_out);
                //fprintf(stderr, "  calling freedv_comptx_700d()\n");
                f->modem_frame_count_tx = 0;
            }
            /* output n_nom_modem_samples at a time from modulated buffer */
            for(i=0; i<f->n_nat_modem_samples; i++) {
                mod_out[i] = 
		    f->mod_out[f->modem_frame_count_tx * f->n_nat_modem_samples+i];
            }
	}
    }
    
#ifdef __LPCNET__
    /* special treatment due to interleaver */
    
    if (FDV_MODE_ACTIVE( FREEDV_MODE_2020, f->mode)) {
        int bits_per_codec_frame = lpcnet_bits_per_frame(f->lpcnet);
        int data_bits_per_frame = f->ldpc->data_bits_per_frame;
	int codec_frames = data_bits_per_frame / bits_per_codec_frame;
        
        //fprintf(stderr, "modem_frame_count_tx: %d dec_frames: %d bytes offset: %d\n",
        //        f->modem_frame_count_tx, codec_frames, (f->modem_frame_count_tx*codec_frames)*bytes_per_codec_frame);
       
        /* buffer up bits until we get enough encoded bits for interleaver */
        
        for (j=0; j<codec_frames; j++) {
            lpcnet_enc(f->lpcnet, speech_in, (char*)f->packed_codec_bits + (f->modem_frame_count_tx*codec_frames+j)*bits_per_codec_frame);
            speech_in += lpcnet_samples_per_frame(f->lpcnet);
        }

	/* Only use extra local buffer if needed for interleave > 1 */
	if (f->interleave_frames == 1) {
            freedv_comptx_2020(f, mod_out);
	} else {
            /* call modulate function when we have enough frames to run interleaver */
            assert((f->modem_frame_count_tx >= 0) && 
	    		(f->modem_frame_count_tx < f->interleave_frames));
            f->modem_frame_count_tx++;
            if (f->modem_frame_count_tx == f->interleave_frames) {
                freedv_comptx_2020(f, f->mod_out);
                //fprintf(stderr, "  calling freedv_comptx_700d()\n");
                f->modem_frame_count_tx = 0;
            }
            /* output n_nom_modem_samples at a time from modulated buffer */
            for(i=0; i<f->n_nat_modem_samples; i++) {
                mod_out[i] = 
		    f->mod_out[f->modem_frame_count_tx * f->n_nat_modem_samples+i];
            }
	}
    }
#endif
    
    /* 2400 A and B are handled by the real-mode TX */
    if(FDV_MODE_ACTIVE( FREEDV_MODE_2400A, f->mode) || FDV_MODE_ACTIVE( FREEDV_MODE_2400B, f->mode)){
    	codec2_encode(f->codec2, f->packed_codec_bits, speech_in);
        freedv_comptx_fsk_voice(f,mod_out);
    }
}

void freedv_codectx(struct freedv *f, short mod_out[], unsigned char *packed_codec_bits) {
    assert(f != NULL);
    COMP tx_fdm[f->n_nom_modem_samples];
    int bits_per_codec_frame;
    int bytes_per_codec_frame;
    int codec_frames;
    int  i;
    bits_per_codec_frame = codec2_bits_per_frame(f->codec2);
    bytes_per_codec_frame = (bits_per_codec_frame + 7) / 8;
    codec_frames = f->n_codec_bits / bits_per_codec_frame;

    memcpy(f->packed_codec_bits, packed_codec_bits, bytes_per_codec_frame * codec_frames);
    
    switch(f->mode) {
        case FREEDV_MODE_1600:
            freedv_comptx_fdmdv_1600(f, tx_fdm);
            break;
        case FREEDV_MODE_700:
        case FREEDV_MODE_700B:
        case FREEDV_MODE_700C:
            freedv_comptx_700(f, tx_fdm);
            break;
        case FREEDV_MODE_700D: {
            /* special treatment due to interleaver */
            int data_bits_per_frame = f->ldpc->data_bits_per_frame;
	    int codec_frames = data_bits_per_frame / bits_per_codec_frame;
	    int j;

            /* buffer up bits until we get enough encoded bits for interleaver */
        
            for (j=0; j<codec_frames; j++) {
                memcpy(f->packed_codec_bits_tx + (f->modem_frame_count_tx*codec_frames+j)*bytes_per_codec_frame, packed_codec_bits, bytes_per_codec_frame);
	        packed_codec_bits += bytes_per_codec_frame;
            }

            /* call modulate function when we have enough frames to run interleaver */

            assert((f->modem_frame_count_tx >= 0) && (f->modem_frame_count_tx < f->interleave_frames));
            f->modem_frame_count_tx++;
            if (f->modem_frame_count_tx == f->interleave_frames) {
                freedv_comptx_700d(f, f->mod_out);
                f->modem_frame_count_tx = 0;
            }

            /* output n_nom_modem_samples at a time from modulated buffer */
            for(i=0; i<f->n_nat_modem_samples; i++) {
                mod_out[i] = f->mod_out[f->modem_frame_count_tx*f->n_nat_modem_samples+i].real;
            }

	    return; /* output is already real */
	}
        case FREEDV_MODE_2400A:
        case FREEDV_MODE_2400B:
        case FREEDV_MODE_800XA:
            freedv_tx_fsk_voice(f, mod_out);
            return; /* output is already real */
    }
    /* convert complex to real */
    for(i=0; i<f->n_nom_modem_samples; i++)
        mod_out[i] = tx_fdm[i].real;
}

void freedv_datatx  (struct freedv *f, short mod_out[]){
    assert(f != NULL);
    if (FDV_MODE_ACTIVE( FREEDV_MODE_2400A, f->mode) || FDV_MODE_ACTIVE( FREEDV_MODE_2400B, f->mode) || FDV_MODE_ACTIVE( FREEDV_MODE_800XA, f->mode)) {
            freedv_tx_fsk_data(f, mod_out);
    }
}

int  freedv_data_ntxframes (struct freedv *f){
    assert(f != NULL);
    if (FDV_MODE_ACTIVE( FREEDV_MODE_2400A, f->mode) || FDV_MODE_ACTIVE( FREEDV_MODE_2400B, f->mode)) {
        if (f->deframer->fdc)
            return freedv_data_get_n_tx_frames(f->deframer->fdc, 8);
    } else if (FDV_MODE_ACTIVE( FREEDV_MODE_800XA, f->mode)) {
        if (f->deframer->fdc)
            return freedv_data_get_n_tx_frames(f->deframer->fdc, 6);
    }
    return 0;
}

int freedv_nin(struct freedv *f) {
    if ((FDV_MODE_ACTIVE( FREEDV_MODE_700, f->mode)) || (FDV_MODE_ACTIVE( FREEDV_MODE_700B, f->mode)) || (FDV_MODE_ACTIVE( FREEDV_MODE_700C, f->mode)))
        // For mode 700, the input rate is 8000 sps, but the modem rate is 7500 sps
        // For mode 700, we request a larger number of Rx samples that will be decimated to f->nin samples
        return (16 * f->nin + f->ptFilter8000to7500->decim_index) / 15;
    else
        return f->nin;
}

/*---------------------------------------------------------------------------*\

  FUNCTION....: freedv_rx
  AUTHOR......: David Rowe
  DATE CREATED: 3 August 2014

  Takes a frame of samples from the radio receiver, demodulates and
  decodes them, producing a frame of decoded speech samples.  See
  freedv_rx.c for an example.

  demod_in[] is a block of received samples sampled at 8000 Hz.
  To account for difference in the transmit and receive sample clock
  frequencies, the number of demod_in[] samples is time varying. You
  MUST call freedv_nin() BEFORE each call to freedv_rx() and pass
  exactly that many samples to this function.

  To help set your buffer sizes, The maximum value of freedv_nin() is
  freedv_get_n_max_modem_samples().

  freedv_rx() returns the number of output speech samples available in
  speech_out[], which is sampled at 8000 Hz.  You should ALWAYS check
  the return value of freedv_rx(), and read EXACTLY that number of
  speech samples from speech_out[].

  1600 and 700D mode: When out of sync, the number of output speech
  samples returned will be freedv_nin(). When in sync to a valid
  FreeDV 1600 signal, the number of output speech samples will
  alternate between freedv_get_n_speech_samples() and 0.

  700 .. 700C modes: The number of output speech samples returned will
  always be freedv_get_n_speech_samples(), regardless of sync.

  The peak level of demod_in[] is not critical, as the demod works
  well over a wide range of amplitude scaling.  However avoid clipping
  (overload, or samples pinned to +/- 32767).  speech_out[] will peak
  at just less than +/-32767.

  When out of sync, this function echoes the demod_in[] samples to
  speech_out[].  This allows the user to listen to the channel, which
  is useful for tuning FreeDV signals or reception of non-FreeDV
  signals.  Setting the squelch with freedv_set_squelch_en(1) will
  return zero-valued samples instead.

\*---------------------------------------------------------------------------*/


// short version

int freedv_rx(struct freedv *f, short speech_out[], short demod_in[]) {
    assert(f != NULL);
    int i;
    int nin = freedv_nin(f);

    assert(nin <= f->n_max_modem_samples);
       
    /* FSK RX happens in real floats, so convert to those and call their demod here */
    if( (FDV_MODE_ACTIVE( FREEDV_MODE_2400A, f->mode)) || (FDV_MODE_ACTIVE( FREEDV_MODE_2400B, f->mode)) || (FDV_MODE_ACTIVE( FREEDV_MODE_800XA, f->mode)) ){
        float rx_float[f->n_max_modem_samples];
        for(i=0; i<nin; i++) {
            rx_float[i] = ((float)demod_in[i]);
        }
        return freedv_floatrx(f,speech_out,rx_float);
    }
    
    if ( (FDV_MODE_ACTIVE( FREEDV_MODE_1600, f->mode)) || (FDV_MODE_ACTIVE( FREEDV_MODE_700, f->mode)) || 
         (FDV_MODE_ACTIVE( FREEDV_MODE_700B, f->mode)) || (FDV_MODE_ACTIVE( FREEDV_MODE_700C, f->mode))
         || (FDV_MODE_ACTIVE( FREEDV_MODE_2020, f->mode))) {

        float gain = 1.0;
        
        assert(nin <= f->n_max_modem_samples);
        COMP rx_fdm[f->n_max_modem_samples];
        for(i=0; i<nin; i++) {
            rx_fdm[i].real = gain*(float)demod_in[i];
            rx_fdm[i].imag = 0.0;
        }
        return freedv_comprx(f, speech_out, rx_fdm);
    }

    /* special low memory version for 700D, to help with stm32 port */
    if (FDV_MODE_ACTIVE( FREEDV_MODE_700D, f->mode)) {
        float gain = 2.0; /* keep levels the same as Octave simulations and C unit tests for real signals */
        return freedv_shortrx(f, speech_out, demod_in, gain);
    }
    
    return 0; /* should never get here */
}


// float input samples version
int freedv_comprx_fsk(struct freedv *f, COMP demod_in[], int *valid) {
    /* Varicode and protocol bits */
    uint8_t vc_bits[2];
    uint8_t proto_bits[3];
    short vc_bit;
    int i;
    int n_ascii;
    char ascii_out;

    if(FDV_MODE_ACTIVE( FREEDV_MODE_2400A, f->mode) || FDV_MODE_ACTIVE( FREEDV_MODE_800XA, f->mode)){        
	fsk_demod(f->fsk,(uint8_t*)f->tx_bits,demod_in);
        f->nin = fsk_nin(f->fsk);
        float EbNodB = f->fsk->stats->snr_est;           /* fsk demod actually estimates Eb/No     */
        f->snr_est = EbNodB + 10.0*log10f(800.0/3000.0); /* so convert to SNR Rb=800, noise B=3000 */
        //fprintf(stderr," %f %f\n", EbNodB, f->snr_est);
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
    
    if(fvhff_deframe_bits(f->deframer,f->packed_codec_bits,proto_bits,vc_bits,(uint8_t*)f->tx_bits)){
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
        *valid = 1;

        /* squelch if if sync but SNR too low */
        if (f->squelch_en && (f->snr_est < f->snr_squelch_thresh)) {
            *valid = 0;
        }
    } else {
        /* squelch if out of sync, or echo input of squelch off */
        if (f->squelch_en) 
            *valid = 0;
        else
            *valid = -1;
    }
    f->sync = f->deframer->state;
    f->stats.sync = f->deframer->state;

    return f->n_speech_samples;
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

// complex input samples version

static int freedv_comprx_fdmdv_1600(struct freedv *f, COMP demod_in[], int *valid) {
    int                 bits_per_codec_frame, bytes_per_codec_frame, bits_per_fdmdv_frame;
    int                 i, j, bit, byte, nin_prev, nout;
    int                 recd_codeword, codeword1, data_flag_index, n_ascii;
    short               abit[1];
    char                ascii_out;
    int                 reliable_sync_bit;

    bits_per_codec_frame  = codec2_bits_per_frame(f->codec2);
    bytes_per_codec_frame = (bits_per_codec_frame + 7) / 8;
    nout = f->n_speech_samples;

    COMP ademod_in[f->nin];
    for(i=0; i<f->nin; i++)
        ademod_in[i] = fcmult(1.0/FDMDV_SCALE, demod_in[i]);

    bits_per_fdmdv_frame  = fdmdv_bits_per_frame(f->fdmdv);

    nin_prev = f->nin;
    fdmdv_demod(f->fdmdv, f->fdmdv_bits, &reliable_sync_bit, ademod_in, &f->nin);
    fdmdv_get_demod_stats(f->fdmdv, &f->stats);
    f->sync = f->fdmdv->sync;
    f->snr_est = f->stats.snr_est;

    if (reliable_sync_bit == 1) {
        f->evenframe = 1;
    }

    if (f->stats.sync) {
        if (f->evenframe == 0) {
            memcpy(f->rx_bits, f->fdmdv_bits, bits_per_fdmdv_frame*sizeof(int));
            nout = 0;
	    *valid = 0;
        }
        else {
            memcpy(&f->rx_bits[bits_per_fdmdv_frame], f->fdmdv_bits, bits_per_fdmdv_frame*sizeof(int));

            if (f->test_frames == 0) {
                recd_codeword = 0;
                for(i=0; i<8; i++) {
                    recd_codeword <<= 1;
                    recd_codeword |= (f->rx_bits[i] & 0x1);
                }
                for(i=11; i<15; i++) {
                    recd_codeword <<= 1;
                    recd_codeword |= (f->rx_bits[i] & 0x1);
                }
                for(i=bits_per_codec_frame; i<bits_per_codec_frame+11; i++) {
                    recd_codeword <<= 1;
                    recd_codeword |= (f->rx_bits[i] & 0x1);
                }
                codeword1 = golay23_decode(recd_codeword);
                f->total_bit_errors += golay23_count_errors(recd_codeword, codeword1);
                f->total_bits       += 23;

                //codeword1 = recd_codeword;
                //fprintf(stderr, "received codeword1: 0x%x  decoded codeword1: 0x%x\n", recd_codeword, codeword1);

                for(i=0; i<bits_per_codec_frame; i++)
                    f->codec_bits[i] = f->rx_bits[i];

                for(i=0; i<8; i++) {
                    f->codec_bits[i] = (codeword1 >> (22-i)) & 0x1;
                }
                for(i=8,j=11; i<12; i++,j++) {
                    f->codec_bits[j] = (codeword1 >> (22-i)) & 0x1;
                }

                // extract txt msg data bit ------------------------------------------------------------

                data_flag_index = codec2_get_spare_bit_index(f->codec2);
                abit[0] = f->codec_bits[data_flag_index];

                n_ascii = varicode_decode(&f->varicode_dec_states, &ascii_out, abit, 1, 1);
                if (n_ascii && (f->freedv_put_next_rx_char != NULL)) {
                    (*f->freedv_put_next_rx_char)(f->callback_state, ascii_out);
                }

                // reconstruct missing bit we steal for data bit and decode speech

                codec2_rebuild_spare_bit(f->codec2, f->codec_bits);

                // pack bits, MSB received first

                bit  = 7;
                byte = 0;
                memset(f->packed_codec_bits, 0,  bytes_per_codec_frame);
                for(i=0; i<bits_per_codec_frame; i++) {
                    f->packed_codec_bits[byte] |= (f->codec_bits[i] << bit);
                    bit--;
                    if(bit < 0) {
                        bit = 7;
                        byte++;
                    }
                }
                *valid = 1;
            }
            else {
                int   test_frame_sync, bit_errors, ntest_bits, k;
                short error_pattern[fdmdv_error_pattern_size(f->fdmdv)];

                for(k=0; k<2; k++) {
                    /* test frames, so lets sync up to the test frames and count any errors */

                    fdmdv_put_test_bits(f->fdmdv, &test_frame_sync, error_pattern, &bit_errors, &ntest_bits, &f->rx_bits[k*bits_per_fdmdv_frame]);

                    if (test_frame_sync == 1) {
                        f->test_frame_sync_state = 1;
                        f->test_frame_count = 0;
                    }

                    if (f->test_frame_sync_state) {
                        if (f->test_frame_count == 0) {
                            f->total_bit_errors += bit_errors;
                            f->total_bits += ntest_bits;
                            if (f->freedv_put_error_pattern != NULL) {
                                (*f->freedv_put_error_pattern)(f->error_pattern_callback_state, error_pattern, fdmdv_error_pattern_size(f->fdmdv));
                            }
                        }
                        f->test_frame_count++;
                        if (f->test_frame_count == 4)
                            f->test_frame_count = 0;
                    }

                    //fprintf(stderr, "test_frame_sync: %d test_frame_sync_state: %d bit_errors: %d ntest_bits: %d\n",
                    //        test_frame_sync, f->test_frame_sync_state, bit_errors, ntest_bits);
                }
            }


            /* squelch if beneath SNR threshold or test frames enabled */

            if ((f->squelch_en && (f->stats.snr_est < f->snr_squelch_thresh)) || f->test_frames) {
                //fprintf(stderr,"squelch %f %f !\n", f->stats.snr_est, f->snr_squelch_thresh);
                *valid = 0;
            }

            nout = f->n_speech_samples;

        }

        /* note this freewheels if reliable sync dissapears on bad channels */

        if (f->evenframe)
            f->evenframe = 0;
        else
            f->evenframe = 1;
        //fprintf(stderr,"%d\n",  f->evenframe);

    } /* if (sync) .... */
    else {
        /* if not in sync pass through analog samples */
        /* this lets us "hear" whats going on, e.g. during tuning */

        //fprintf(stderr, "out of sync\n");

        if (f->squelch_en == 0) {
	    *valid = -1;
        }
        else {
	    *valid = 0;
        }
        //fprintf(stderr, "%d %d %d\n", nin_prev, speech_out[0], speech_out[nin_prev-1]);
        nout = nin_prev;
    }
    return nout;
}

static int freedv_comprx_700(struct freedv *f, COMP demod_in_8kHz[], int *valid) {
    int                 bits_per_codec_frame, bytes_per_codec_frame;
    int                 i, j, bit, byte, nout, k;
    int                 data_flag_index, n_ascii, nspare;
    short               abit[1];
    char                ascii_out;
    float rx_bits[COHPSK_BITS_PER_FRAME]; /* soft decn rx bits */
    int   sync;
    int   frames;

    bits_per_codec_frame  = codec2_bits_per_frame(f->codec2);
    bytes_per_codec_frame = (bits_per_codec_frame + 7) / 8;
    frames = f->n_codec_bits / bits_per_codec_frame;
    nout = f->n_speech_samples;

    // echo samples back out as default (say if sync not found)
    *valid = -1;

    // quisk_cfInterpDecim() modifies input data so lets make a copy just in case there
    // is no sync and we need to echo inout to output

    COMP demod_in[freedv_nin(f)];
    for(i=0; i<freedv_nin(f); i++)
        demod_in[i] = demod_in_8kHz[i];

    i = quisk_cfInterpDecim((complex float *)demod_in, freedv_nin(f), f->ptFilter8000to7500, 15, 16);
    //if (i != f->nin)
    //    printf("freedv_comprx decimation: input %d output %d\n", freedv_nin(f), i);

    for(i=0; i<f->nin; i++)
        demod_in[i] = fcmult(1.0/FDMDV_SCALE, demod_in[i]);
    
    cohpsk_demod(f->cohpsk, rx_bits, &sync, demod_in, &f->nin);

    f->sync = sync;
    cohpsk_get_demod_stats(f->cohpsk, &f->stats);
    f->snr_est = f->stats.snr_est;

    memset(f->packed_codec_bits, 0, bytes_per_codec_frame * frames);

    if (sync) {

        if (f->test_frames == 0) {
            data_flag_index = codec2_get_spare_bit_index(f->codec2);

            /* optional smoothing of codec symbols */

            if (f->smooth_symbols) {

                for(i=0; i<bits_per_codec_frame; i++) {
                    rx_bits[i] += rx_bits[i+bits_per_codec_frame];
                    rx_bits[i+bits_per_codec_frame] = rx_bits[i];
                }
            }

            byte = 0;
            for (j=0; j<COHPSK_BITS_PER_FRAME; j+=bits_per_codec_frame) {

                /* extract txt msg data bit(s) */

                switch(f->mode) {
                case FREEDV_MODE_700:
                    nspare = 2;
                    break;
                case FREEDV_MODE_700B:
                    nspare = 1; // Just one spare bit for FREEDV_MODE_700B
                    break;
                case FREEDV_MODE_700C:
                    nspare = 0; // and no spare bits for 700C atm
                    break;
                default:
                    assert(0);
                }

                for(k=0; k<nspare; k++)  {
                    abit[0] = rx_bits[data_flag_index+j+k] < 0.0;

                    n_ascii = varicode_decode(&f->varicode_dec_states, &ascii_out, abit, 1, 1);
                    if (n_ascii && (f->freedv_put_next_rx_char != NULL)) {
                        (*f->freedv_put_next_rx_char)(f->callback_state, ascii_out);
                    }
                }

                /* pack bits, MSB received first */

                bit = 7;
                for(i=0; i<bits_per_codec_frame; i++) {
                    f->packed_codec_bits[byte] |= ((rx_bits[j+i] < 0.0) << bit);
                    bit--;
                    if (bit < 0) {
                        bit = 7;
                        byte++;
                    }
                }
		if (bit != 7)
		    byte++;

                if (f->squelch_en && (f->stats.snr_est < f->snr_squelch_thresh)) {
		   *valid = 0;
                }
		*valid = 1;
            }
            nout = f->n_speech_samples;
        }
        else {
            //fprintf(stderr, " freedv_api:  f->test_frames_diversity: %d\n", f->test_frames_diversity);

            if (f->test_frames_diversity) {
                /* normal operation - error pattern on frame after diveristy combination */
                short error_pattern[COHPSK_BITS_PER_FRAME];
                int   bit_errors;

                /* test data, lets see if we can sync to the test data sequence */

                char rx_bits_char[COHPSK_BITS_PER_FRAME];
                for(i=0; i<COHPSK_BITS_PER_FRAME; i++)
                    rx_bits_char[i] = rx_bits[i] < 0.0;
                cohpsk_put_test_bits(f->cohpsk, &f->test_frame_sync_state, error_pattern, &bit_errors, rx_bits_char, 0);
                if (f->test_frame_sync_state) {
                    f->total_bit_errors += bit_errors;
                    f->total_bits       += COHPSK_BITS_PER_FRAME;
                    if (f->freedv_put_error_pattern != NULL) {
                        (*f->freedv_put_error_pattern)(f->error_pattern_callback_state, error_pattern, COHPSK_BITS_PER_FRAME);
                    }
                }
            } 
            else {
                /* calculate error pattern on uncombined carriers - test mode to spot any carrier specific issues like
                   tx passband filtering */

                short error_pattern[2*COHPSK_BITS_PER_FRAME];
                char  rx_bits_char[COHPSK_BITS_PER_FRAME];
                int   bit_errors_lower, bit_errors_upper;

                /* lower group of carriers */

                float *rx_bits_lower = cohpsk_get_rx_bits_lower(f->cohpsk);
                for(i=0; i<COHPSK_BITS_PER_FRAME; i++) {
                    rx_bits_char[i] = rx_bits_lower[i] < 0.0;
                }
                cohpsk_put_test_bits(f->cohpsk, &f->test_frame_sync_state, error_pattern, &bit_errors_lower, rx_bits_char, 0);

                /* upper group of carriers */

                float *rx_bits_upper = cohpsk_get_rx_bits_upper(f->cohpsk);
                for(i=0; i<COHPSK_BITS_PER_FRAME; i++) {
                    rx_bits_char[i] = rx_bits_upper[i] < 0.0;
                }
                cohpsk_put_test_bits(f->cohpsk, &f->test_frame_sync_state_upper, &error_pattern[COHPSK_BITS_PER_FRAME], &bit_errors_upper, rx_bits_char, 1);
                //                fprintf(stderr, " freedv_api:  f->test_frame_sync_state: %d f->test_frame_sync_state_upper: %d\n", 
                //        f->test_frame_sync_state, f->test_frame_sync_state_upper);

                /* combine total errors and call callback */

                if (f->test_frame_sync_state && f->test_frame_sync_state_upper) {
                    f->total_bit_errors += bit_errors_lower + bit_errors_upper;
                    f->total_bits       += 2*COHPSK_BITS_PER_FRAME;
                    if (f->freedv_put_error_pattern != NULL) {
                        (*f->freedv_put_error_pattern)(f->error_pattern_callback_state, error_pattern, 2*COHPSK_BITS_PER_FRAME);
                    }
                }

            }
            
	    *valid = 0;
            nout = f->n_speech_samples;
        }

    }

    /* no valid FreeDV signal - squelch output */

    if (sync == 0) {
        nout = freedv_nin(f);
        if (f->squelch_en) {
	    *valid = 0;
        }
    }
    return nout;
}

/*
  700D demod function that can support complex (float) or real (short)
  samples.  The real short samples are useful for low memory overhead,
  such at the SM1000.
*/

static int freedv_comp_short_rx_700d(struct freedv *f, void *demod_in_8kHz, int demod_in_is_short, float gain, int *valid) {
    int   bits_per_codec_frame, bytes_per_codec_frame;
    int   i, j, bit, byte, nout, k;
    int   n_ascii;
    char  ascii_out;
    int   frames;
    struct OFDM *ofdm = f->ofdm;
    struct LDPC *ldpc = f->ldpc;
    
    int    data_bits_per_frame = ldpc->data_bits_per_frame;
    int    coded_bits_per_frame = ldpc->coded_bits_per_frame;
    int    coded_syms_per_frame = ldpc->coded_syms_per_frame;
    int    interleave_frames = f->interleave_frames;
    COMP  *codeword_symbols = f->codeword_symbols;
    float *codeword_amps = f->codeword_amps;
    int    rx_bits[ofdm_bitsperframe];
    short txt_bits[ofdm_ntxtbits];
    COMP  payload_syms[coded_syms_per_frame];
    float payload_amps[coded_syms_per_frame];
   
    bits_per_codec_frame  = codec2_bits_per_frame(f->codec2);
    bytes_per_codec_frame = (bits_per_codec_frame + 7) / 8;
    frames = f->n_codec_bits / bits_per_codec_frame;

    // pass through is too noisey ....
    //nout = f->n_speech_samples;
    nout = 0;
    
    int Nerrs_raw = 0;
    int Nerrs_coded = 0;
    int iter = 0;
    int parityCheckCount = 0;
    uint8_t rx_uw[ofdm_nuwbits];

    float new_gain = gain / OFDM_AMP_SCALE;
    
    assert((demod_in_is_short == 0) || (demod_in_is_short == 1));

    /* echo samples back out as default (say if sync not found) */
    
    *valid = 1;
    f->sync = f->stats.sync = 0;
    
    /* TODO estimate this properly from signal */
    
    float EsNo = 3.0;
    
    /* looking for modem sync */
    
    if (ofdm->sync_state == search) {
        if (demod_in_is_short)
            ofdm_sync_search_shorts(f->ofdm, (short*)demod_in_8kHz, new_gain);
        else
            ofdm_sync_search(f->ofdm, (COMP*)demod_in_8kHz);
    }

     /* OK modem is in sync */
    
    if ((ofdm->sync_state == synced) || (ofdm->sync_state == trial)) {
        if (demod_in_is_short)
            ofdm_demod_shorts(ofdm, rx_bits, (short*)demod_in_8kHz, new_gain);
        else
            ofdm_demod(ofdm, rx_bits, (COMP*)demod_in_8kHz);
            
        ofdm_disassemble_modem_frame(ofdm, rx_uw, payload_syms, payload_amps, txt_bits);

        f->sync = 1;
        ofdm_get_demod_stats(f->ofdm, &f->stats);
        f->snr_est = f->stats.snr_est;

        assert((ofdm_nuwbits+ofdm_ntxtbits+coded_bits_per_frame) == ofdm_bitsperframe);

        /* now we need to buffer for de-interleaving -------------------------------------*/
                
        /* shift interleaved symbol buffers to make room for new symbols */
                
        for(i=0, j=coded_syms_per_frame; j<interleave_frames*coded_syms_per_frame; i++,j++) {
            codeword_symbols[i] = codeword_symbols[j];
            codeword_amps[i] = codeword_amps[j];
        }

        /* newest symbols at end of buffer (uses final i from last loop), note we 
           change COMP formats from what modem uses internally */
                
        for(i=(interleave_frames-1)*coded_syms_per_frame,j=0; i<interleave_frames*coded_syms_per_frame; i++,j++) {
            codeword_symbols[i] = payload_syms[j];
            codeword_amps[i]    = payload_amps[j];
         }
               
        /* run de-interleaver */
                
        COMP  codeword_symbols_de[interleave_frames*coded_syms_per_frame];
        float codeword_amps_de[interleave_frames*coded_syms_per_frame];
        gp_deinterleave_comp (codeword_symbols_de, codeword_symbols, interleave_frames*coded_syms_per_frame);
        gp_deinterleave_float(codeword_amps_de   , codeword_amps   , interleave_frames*coded_syms_per_frame);

        float llr[coded_bits_per_frame];
        uint8_t out_char[coded_bits_per_frame];

        interleaver_sync_state_machine(ofdm, ldpc, ofdm_config, codeword_symbols_de, codeword_amps_de, EsNo,
                                       interleave_frames, &iter, &parityCheckCount, &Nerrs_coded);
                                         
        if ((ofdm->sync_state_interleaver == synced) && (ofdm->frame_count_interleaver == interleave_frames)) {
            ofdm->frame_count_interleaver = 0;

            if (f->test_frames) {
                int tmp[interleave_frames];
                Nerrs_raw = count_uncoded_errors(ldpc, ofdm_config, tmp, interleave_frames, codeword_symbols_de);
                f->total_bit_errors += Nerrs_raw;
                f->total_bits       += ofdm_bitsperframe*interleave_frames;
            }

            memset(f->packed_codec_bits, 0, bytes_per_codec_frame * frames);
            byte = 0; f->modem_frame_count_rx = 0;
            
            for (j=0; j<interleave_frames; j++) {
                symbols_to_llrs(llr, &codeword_symbols_de[j*coded_syms_per_frame],
                                &codeword_amps_de[j*coded_syms_per_frame],
                                EsNo, ofdm->mean_amp, coded_syms_per_frame);               
                iter = run_ldpc_decoder(ldpc, out_char, llr, &parityCheckCount);

                if (f->test_frames) {
                    uint8_t payload_data_bits[data_bits_per_frame];
                    ofdm_generate_payload_data_bits(payload_data_bits, data_bits_per_frame);
                    Nerrs_coded = count_errors(payload_data_bits, out_char, data_bits_per_frame);
                    f->total_bit_errors_coded += Nerrs_coded;
                    f->total_bits_coded       += data_bits_per_frame;
                } else {

                    /* a frame of valid Codec 2 bits, pack into Codec 2 frame  */

                    for (i=0; i<data_bits_per_frame; i+=bits_per_codec_frame) {

                        /* pack bits, MSB received first */

                        bit = 7;
                        for(k=0; k<bits_per_codec_frame; k++) {
                            f->packed_codec_bits[byte] |= (out_char[i+k] << bit);
                            bit--;
                            if (bit < 0) {
                                bit = 7;
                                byte++;
                            }
                        }
                        if (bit != 7)
                            byte++;
                    }
                    
                }
            } /* for interleave frames ... */

            /* make sure we don't overrun packed byte array */

            assert(byte <= f->nbyte_packed_codec_bits);
                   
            nout = f->n_speech_samples;                  

            if (f->squelch_en && (f->stats.snr_est < f->snr_squelch_thresh)) {
                *valid = 0;
            }
            
        } /* if interleaver synced ..... */

        /* If modem is synced we can decode txt bits */
        
        for(k=0; k<ofdm_ntxtbits; k++)  { 
            //fprintf(stderr, "txt_bits[%d] = %d\n", k, rx_bits[i]);
            n_ascii = varicode_decode(&f->varicode_dec_states, &ascii_out, &txt_bits[k], 1, 1);
            if (n_ascii && (f->freedv_put_next_rx_char != NULL)) {
                (*f->freedv_put_next_rx_char)(f->callback_state, ascii_out);
            }
        }

        /* estimate uncoded BER from UW.  Coded bit errors could
           probably be estimated as half of all failed LDPC parity
           checks */

        for(i=0; i<ofdm_nuwbits; i++) {         
            if (rx_uw[i] != ofdm->tx_uw[i]) {
                f->total_bit_errors++;
            }
        }
        f->total_bits += ofdm_nuwbits;          

    } /* if modem synced .... */ else {
        *valid = -1;
    }

    /* iterate state machine and update nin for next call */
    
    f->nin = ofdm_get_nin(ofdm);
    ofdm_sync_state_machine(ofdm, rx_uw);

    /* check if OFDM modem has sync */    
    bool sync = ((ofdm->sync_state == synced) || (ofdm->sync_state == trial));

    /* if LDPC code has too many errors it may be the end of an over or a very deep fade */
    bool ldpc_decode_ok = parityCheckCount > 0.8*ldpc->NumberParityBits;

    /* squelch if out of sync or too many LDPC decode errors */
    if ((sync == false) || (ldpc_decode_ok == false)) {
        if (f->squelch_en == true) {
 	    *valid = 0;
        }
    }
        
    if ((f->verbose && (ofdm->last_sync_state == search)) || (f->verbose == 2)) {
        fprintf(stderr, "%3d nin: %4d st: %-6s euw: %2d %1d f: %5.1f phbw: %d ist: %-6s %2d eraw: %3d ecdd: %3d iter: %3d pcc: %3d vld: %d, nout: %4d\n",
                f->frames++, ofdm->nin, statemode[ofdm->last_sync_state], ofdm->uw_errors, ofdm->sync_counter, 
		(double)ofdm->foff_est_hz, ofdm->phase_est_bandwidth,
                statemode[ofdm->last_sync_state_interleaver], ofdm->frame_count_interleaver,
                Nerrs_raw, Nerrs_coded, iter, parityCheckCount, *valid, nout);
    }
    
    return nout;
}

#ifdef __LPCNET__
static int freedv_comprx_2020(struct freedv *f, COMP demod_in[], int *valid) {
    int   i, j, nout, k;
    int   n_ascii;
    char  ascii_out;
    struct OFDM *ofdm = f->ofdm;
    struct LDPC *ldpc = f->ldpc;
    
    int    data_bits_per_frame = ldpc->data_bits_per_frame;
    int    coded_bits_per_frame = ldpc->coded_bits_per_frame;
    int    coded_syms_per_frame = ldpc->coded_syms_per_frame;
    int    interleave_frames = f->interleave_frames;
    COMP  *codeword_symbols = f->codeword_symbols;
    float *codeword_amps = f->codeword_amps;
    int    rx_bits[ofdm_bitsperframe];
    short txt_bits[ofdm_ntxtbits];
    COMP  payload_syms[coded_syms_per_frame];
    float payload_amps[coded_syms_per_frame];
   
    // pass through is too noisey ....
    //nout = f->n_speech_samples;
    nout = 0;
    
    int Nerrs_raw = 0;
    int Nerrs_coded = 0;
    int iter = 0;
    int parityCheckCount = 0;
    uint8_t rx_uw[ofdm_nuwbits];

    /* echo samples back out as default (say if sync not found) */
    
    *valid = 1;
    f->sync = f->stats.sync = 0;
    
    /* TODO estimate this properly from signal */

    // TODO: should be higher for 2020?
    float EsNo = 3.0;
    
    /* looking for modem sync */
    
    if (ofdm->sync_state == search) {
        ofdm_sync_search(f->ofdm, demod_in);
    }

    /* OK modem is in sync */
    
    if ((ofdm->sync_state == synced) || (ofdm->sync_state == trial)) {
        ofdm_demod(ofdm, rx_bits, demod_in);
        ofdm_disassemble_modem_frame(ofdm, rx_uw, payload_syms, payload_amps, txt_bits);

        f->sync = 1;
        ofdm_get_demod_stats(f->ofdm, &f->stats);
        f->snr_est = f->stats.snr_est;
        
        assert((ofdm_nuwbits+ofdm_ntxtbits+coded_bits_per_frame) == ofdm_bitsperframe);

        /* now we need to buffer for de-interleaving -------------------------------------*/
                
        /* shift interleaved symbol buffers to make room for new symbols */
                
        for(i=0, j=coded_syms_per_frame; j<interleave_frames*coded_syms_per_frame; i++,j++) {
            codeword_symbols[i] = codeword_symbols[j];
            codeword_amps[i] = codeword_amps[j];
        }

        /* newest symbols at end of buffer (uses final i from last loop), note we 
           change COMP formats from what modem uses internally */
                
        for(i=(interleave_frames-1)*coded_syms_per_frame,j=0; i<interleave_frames*coded_syms_per_frame; i++,j++) {
            codeword_symbols[i] = payload_syms[j];
            codeword_amps[i]    = payload_amps[j];
         }
               
        /* run de-interleaver */
                
        COMP  codeword_symbols_de[interleave_frames*coded_syms_per_frame];
        float codeword_amps_de[interleave_frames*coded_syms_per_frame];
        gp_deinterleave_comp (codeword_symbols_de, codeword_symbols, interleave_frames*coded_syms_per_frame);
        gp_deinterleave_float(codeword_amps_de   , codeword_amps   , interleave_frames*coded_syms_per_frame);

        float llr[coded_bits_per_frame];
        uint8_t out_char[coded_bits_per_frame];

        interleaver_sync_state_machine(ofdm, ldpc, ofdm_config, codeword_symbols_de, codeword_amps_de, EsNo,
                                       interleave_frames, &iter, &parityCheckCount, &Nerrs_coded);
                                         
        if ((ofdm->sync_state_interleaver == synced) && (ofdm->frame_count_interleaver == interleave_frames)) {
            ofdm->frame_count_interleaver = 0;

            if (f->test_frames) {
                int tmp[interleave_frames];
                Nerrs_raw = count_uncoded_errors(ldpc, ofdm_config, tmp, interleave_frames, codeword_symbols_de);
                f->total_bit_errors += Nerrs_raw;
                f->total_bits       += ofdm_bitsperframe*interleave_frames;
            }

            f->modem_frame_count_rx = 0;
            
            for (j=0; j<interleave_frames; j++) {
                symbols_to_llrs(llr, &codeword_symbols_de[j * coded_syms_per_frame],
                                &codeword_amps_de[j * coded_syms_per_frame],
                                EsNo, ofdm->mean_amp, coded_syms_per_frame);                           
               if (ldpc->data_bits_per_frame == ldpc->ldpc_data_bits_per_frame) {
                    /* all data bits in code word used */
                    iter = run_ldpc_decoder(ldpc, out_char, llr, &parityCheckCount);
                } else {
                    /* some unused data bits, set these to known values to strengthen code */
                    float llr_full_codeword[ldpc->ldpc_coded_bits_per_frame];
                    int unused_data_bits = ldpc->ldpc_data_bits_per_frame - ldpc->data_bits_per_frame;
                                
                    // received data bits
                    for (i = 0; i < ldpc->data_bits_per_frame; i++)
                        llr_full_codeword[i] = llr[i]; 
                    // known bits ... so really likely
                    for (i = ldpc->data_bits_per_frame; i < ldpc->ldpc_data_bits_per_frame; i++)
                        llr_full_codeword[i] = -100.0;
                    // parity bits at end
                    for (i = ldpc->ldpc_data_bits_per_frame; i < ldpc->ldpc_coded_bits_per_frame; i++)
                        llr_full_codeword[i] = llr[i-unused_data_bits]; 
                    iter = run_ldpc_decoder(ldpc, out_char, llr_full_codeword, &parityCheckCount);
                }
                            
                if (f->test_frames) {
                    uint8_t payload_data_bits[data_bits_per_frame];
                    ofdm_generate_payload_data_bits(payload_data_bits, data_bits_per_frame);
                    Nerrs_coded = count_errors(payload_data_bits, out_char, data_bits_per_frame);
                    f->total_bit_errors_coded += Nerrs_coded;
                    f->total_bits_coded       += data_bits_per_frame;
                } else {
                    memcpy(f->packed_codec_bits+j*data_bits_per_frame, out_char, data_bits_per_frame);
                }
            } /* for interleave frames ... */
                   
            nout = f->n_speech_samples;                  

            if (f->squelch_en && (f->stats.snr_est < f->snr_squelch_thresh)) {
                *valid = 0;
            }
            
        } /* if interleaver synced ..... */

        /* If modem is synced we can decode txt bits */
        
        for(k=0; k<ofdm_ntxtbits; k++)  { 
            //fprintf(stderr, "txt_bits[%d] = %d\n", k, rx_bits[i]);
            n_ascii = varicode_decode(&f->varicode_dec_states, &ascii_out, &txt_bits[k], 1, 1);
            if (n_ascii && (f->freedv_put_next_rx_char != NULL)) {
                (*f->freedv_put_next_rx_char)(f->callback_state, ascii_out);
            }
        }

        /* estimate uncoded BER from UW.  Coded bit errors could
           probably be estimated as half of all failed LDPC parity
           checks */

        for(i=0; i<ofdm_nuwbits; i++) {         
            if (rx_uw[i] != ofdm->tx_uw[i]) {
                f->total_bit_errors++;
            }
        }
        f->total_bits += ofdm_nuwbits;          

    } /* if modem synced .... */ else {
        *valid = -1;
    }

    /* iterate state machine and update nin for next call */
    
    f->nin = ofdm_get_nin(ofdm);
    //fprintf(stderr, "nin: %d\n", ofdm_get_nin(ofdm));
    ofdm_sync_state_machine(ofdm, rx_uw);

    if ((f->verbose && (ofdm->last_sync_state == search)) || (f->verbose == 2)) {
        fprintf(stderr, "%3d st: %-6s euw: %2d %1d f: %5.1f pbw: %d ist: %-6s %2d eraw: %3d ecdd: %3d iter: %3d pcc: %3d vld: %d, nout: %4d\n",
                f->frames++, statemode[ofdm->last_sync_state], ofdm->uw_errors, ofdm->sync_counter, 
		(double)ofdm->foff_est_hz, ofdm->phase_est_bandwidth,
                statemode[ofdm->last_sync_state_interleaver], ofdm->frame_count_interleaver,
                Nerrs_raw, Nerrs_coded, iter, parityCheckCount, *valid, nout);
    }
    
    /* check if OFDM modem has sync */    
    bool sync = ((ofdm->sync_state == synced) || (ofdm->sync_state == trial));

    /* if LDPC code has too many errors it may be the end of an over or a very deep fade */
    bool ldpc_decode_ok = parityCheckCount > 0.8*ldpc->NumberParityBits;

    /* squelch if out of sync or too many LDPC decode errors */
    if ((sync == false) || (ldpc_decode_ok == false)) {
         if (f->squelch_en == true) {
 	    *valid = 0;
         }
         //f->snr_est = 0.0;
    }
    
    //fprintf(stderr, "sync: %d valid: %d snr: %3.2f\n", f->sync, *valid, f->snr_est);
    
    return nout;
}
#endif

/* complex input rx API function */
int freedv_comprx(struct freedv *f, short speech_out[], COMP demod_in[]) {
    assert(f != NULL);
    int bits_per_codec_frame, bytes_per_codec_frame;
    int i, nout = 0;
    int valid = 0;
    
    assert(f->nin <= f->n_max_modem_samples);

    if (FDV_MODE_ACTIVE( FREEDV_MODE_1600, f->mode)) {
        nout = freedv_comprx_fdmdv_1600(f, demod_in, &valid);
    }
    if ((FDV_MODE_ACTIVE( FREEDV_MODE_700, f->mode)) || (FDV_MODE_ACTIVE( FREEDV_MODE_700B, f->mode)) || (FDV_MODE_ACTIVE( FREEDV_MODE_700C, f->mode))) {
        nout = freedv_comprx_700(f, demod_in, &valid);
    }

    if( (FDV_MODE_ACTIVE( FREEDV_MODE_2400A, f->mode)) || (FDV_MODE_ACTIVE( FREEDV_MODE_2400B, f->mode)) || (FDV_MODE_ACTIVE( FREEDV_MODE_800XA, f->mode))){
        nout = freedv_comprx_fsk(f, demod_in, &valid);
    }

    if (FDV_MODE_ACTIVE( FREEDV_MODE_700D, f->mode)) {
        nout = freedv_comp_short_rx_700d(f, (void*)demod_in, 0, 1.0, &valid);
    }
    
    if (FDV_MODE_ACTIVE( FREEDV_MODE_2020, f->mode)) {
#ifdef __LPCNET__
        freedv_comprx_2020(f, demod_in, &valid);
#endif
    }
    
    if (valid == 0) {
        /* squelch */
        for (i = 0; i < nout; i++)
            speech_out[i] = 0;
    }
    else if (valid < 0) {
        /* we havent got sync so play undemodulated audio from radio.  This might
           not work for all modes due to nin bouncing about */
        for (i = 0; i < nout; i++)
            speech_out[i] = demod_in[i].real;
    }
    else {
        /* decode audio  -----------------------------------------------*/

        if(FDV_MODE_ACTIVE( FREEDV_MODE_2020, f->mode)) {
#ifdef __LPCNET__
            /* LPCNet decoder */
            int bits_per_codec_frame = lpcnet_bits_per_frame(f->lpcnet);
            int data_bits_per_frame = f->ldpc->data_bits_per_frame;
            int frames = data_bits_per_frame/bits_per_codec_frame;            
        
            /*
              fprintf(stderr, "  modem_frame_count_rx: %d nout: %d frames: %d data_bits_per_frame: %d interleave_frames: %d\n",
                    f->modem_frame_count_rx, nout, frames, data_bits_per_frame, f->interleave_frames);
            */
            if (f->modem_frame_count_rx < f->interleave_frames) {
                nout = f->n_speech_samples;
                for (i = 0; i < frames; i++) {
                    lpcnet_dec(f->lpcnet, (char*)f->packed_codec_bits + (i + frames*f->modem_frame_count_rx)* bits_per_codec_frame, speech_out);
                    speech_out += lpcnet_samples_per_frame(f->lpcnet);
                }
                f->modem_frame_count_rx++;
            }
 #endif          
        }
        else {
            bits_per_codec_frame  = codec2_bits_per_frame(f->codec2);
            bytes_per_codec_frame = (bits_per_codec_frame + 7) / 8;

            if(FDV_MODE_ACTIVE( FREEDV_MODE_700D, f->mode)) {

                /* 700D a bit special due to interleaving */

                if (f->modem_frame_count_rx < f->interleave_frames) {
                    int data_bits_per_frame = f->ldpc->data_bits_per_frame;
                    int frames = data_bits_per_frame/bits_per_codec_frame;            
                    nout = f->n_speech_samples;
                    //fprintf(stderr, "modem_frame_count_rx: %d nout: %d\n", f->modem_frame_count_rx, nout);
                    for (i = 0; i < frames; i++) {
                        codec2_decode(f->codec2, speech_out, f->packed_codec_bits + (i + frames*f->modem_frame_count_rx)* bytes_per_codec_frame);
                        speech_out += codec2_samples_per_frame(f->codec2);
                    }
                    f->modem_frame_count_rx++;
                }
            } else { 
                /* non-interleaved Codec 2 modes */
                
                int frames = f->n_codec_bits / bits_per_codec_frame;
                //fprintf(stderr, "frames: %d\n", frames);
                for (i = 0; i < frames; i++) {
                    codec2_decode(f->codec2, speech_out, f->packed_codec_bits + i * bytes_per_codec_frame);
                    speech_out += codec2_samples_per_frame(f->codec2);
                }
            }
        }
    }

    //fprintf(stderr,"freedv_nin(f): %d nout: %d valid: %d\n", freedv_nin(f), nout, valid);
    return nout;
}

/* memory efficient real short version - just for 700D  */
int freedv_shortrx(struct freedv *f, short speech_out[], short demod_in[], float gain) {
    assert(f != NULL);
    int                 bits_per_codec_frame, bytes_per_codec_frame;
    int                 i, nout = 0;
    int valid = 0;

    // At this stage short interface only supported for 700D, to help
    // memory requirements on stm32
    assert(f->mode == FREEDV_MODE_700D);
    
    assert(f->nin <= f->n_max_modem_samples);

    bits_per_codec_frame  = codec2_bits_per_frame(f->codec2);
    bytes_per_codec_frame = (bits_per_codec_frame + 7) / 8;

    if (FDV_MODE_ACTIVE( FREEDV_MODE_700D, f->mode)) {
        nout = freedv_comp_short_rx_700d(f, (void*)demod_in, 1, gain, &valid);
    }
    
    if (valid == 0) {
        //fprintf(stderr, "squelch nout: %d\n", nout);
        
        /* squelch */
        
        for (i = 0; i < nout; i++)
            speech_out[i] = 0;
    }
    else if (valid < 0) {
        /* we havent got sync so play audio from radio.  This might
           not work for all modes due to nin bouncing about */
        for (i = 0; i < nout; i++)
            speech_out[i] = demod_in[i];
    }
    else {
        /* decoded audio to play */
        
        int data_bits_per_frame = f->ldpc->data_bits_per_frame;
        int frames = data_bits_per_frame/bits_per_codec_frame;
            
        nout = 0;
        if (f->modem_frame_count_rx < f->interleave_frames) {
            nout = f->n_speech_samples;
            //fprintf(stderr, "modem_frame_count_rx: %d nout: %d\n", f->modem_frame_count_rx, nout);
            for (i = 0; i < frames; i++) {
                codec2_decode(f->codec2, speech_out, f->packed_codec_bits + (i + frames*f->modem_frame_count_rx)* bytes_per_codec_frame);
                speech_out += codec2_samples_per_frame(f->codec2);
            }
            f->modem_frame_count_rx++;
        }
    }

    return nout;
}


int freedv_codecrx(struct freedv *f, unsigned char *packed_codec_bits, short demod_in[])
{
    assert(f != NULL);
    int i;
    int nin = freedv_nin(f);
    int valid;
    int ret = 0;
    int bits_per_codec_frame = codec2_bits_per_frame(f->codec2);

    assert(nin <= f->n_max_modem_samples);

    if (FDV_MODE_ACTIVE( FREEDV_MODE_700D, f->mode) == false) {
        COMP rx_fdm[f->n_max_modem_samples];
    
        for(i=0; i<nin; i++) {
            rx_fdm[i].real = (float)demod_in[i];
            rx_fdm[i].imag = 0.0;
        }

        if (FDV_MODE_ACTIVE( FREEDV_MODE_1600, f->mode)) {
            freedv_comprx_fdmdv_1600(f, rx_fdm, &valid);
        }

        if ((FDV_MODE_ACTIVE( FREEDV_MODE_700, f->mode)) || (FDV_MODE_ACTIVE( FREEDV_MODE_700B, f->mode)) || (FDV_MODE_ACTIVE( FREEDV_MODE_700C, f->mode))) {
            freedv_comprx_700(f, rx_fdm, &valid);
        }

        if( FDV_MODE_ACTIVE( FREEDV_MODE_2400A, f->mode) || FDV_MODE_ACTIVE( FREEDV_MODE_2400B, f->mode) || FDV_MODE_ACTIVE( FREEDV_MODE_800XA, f->mode)){
            freedv_comprx_fsk(f, rx_fdm, &valid);
        }
    }

    int bytes_per_codec_frame = (bits_per_codec_frame + 7) / 8;

    if (FDV_MODE_ACTIVE( FREEDV_MODE_700D, f->mode)) {
        freedv_comp_short_rx_700d(f, (void*)demod_in, 1, 1.0, &valid);

        int data_bits_per_frame = f->ldpc->data_bits_per_frame;
        int frames = data_bits_per_frame/bits_per_codec_frame;
            
        if (valid == 1 && f->modem_frame_count_rx < f->interleave_frames) {
             for (i = 0; i < frames; i++) {
                 memcpy(packed_codec_bits, f->packed_codec_bits + (i + frames*f->modem_frame_count_rx)* bytes_per_codec_frame, bytes_per_codec_frame);
                 packed_codec_bits += bytes_per_codec_frame;
                 ret += bytes_per_codec_frame;
             }
             f->modem_frame_count_rx++;
        }
	return ret;
    }

    if (valid == 1) {
        int codec_frames = f->n_codec_bits / bits_per_codec_frame;

        memcpy(packed_codec_bits, f->packed_codec_bits, bytes_per_codec_frame * codec_frames);
	ret = bytes_per_codec_frame * codec_frames;
    }
    
    return ret;
}

/*---------------------------------------------------------------------------*\

  FUNCTION....: freedv_get_version
  AUTHOR......: Jim Ahlstrom
  DATE CREATED: 28 July 2015

  Return the version of the FreeDV API.  This is meant to help API users determine when
  incompatible changes have occurred.

\*---------------------------------------------------------------------------*/

int freedv_get_version(void)
{
    return VERSION;
}

/*---------------------------------------------------------------------------*\

  FUNCTION....: freedv_set_callback_txt
  AUTHOR......: Jim Ahlstrom
  DATE CREATED: 28 July 2015

  Set the callback functions and the callback state pointer that will be used
  for the aux txt channel.  The freedv_callback_rx is a function pointer that
  will be called to return received characters.  The freedv_callback_tx is a
  function pointer that will be called to send transmitted characters.  The callback
  state is a user-defined void pointer that will be passed to the callback functions.
  Any or all can be NULL, and the default is all NULL.
  The function signatures are:
    void receive_char(void *callback_state, char c);
    char transmit_char(void *callback_state);

\*---------------------------------------------------------------------------*/

void freedv_set_callback_txt(struct freedv *f, freedv_callback_rx rx, freedv_callback_tx tx, void *state)
{
    if (FDV_MODE_ACTIVE( FREEDV_MODE_800XA, f->mode) == false) {
        f->freedv_put_next_rx_char = rx;
        f->freedv_get_next_tx_char = tx;
        f->callback_state = state;
    }
}

/*---------------------------------------------------------------------------*\

  FUNCTION....: freedv_set_callback_protocol
  AUTHOR......: Brady OBrien
  DATE CREATED: 21 February 2016

  Set the callback functions and callback pointer that will be used for the
  protocol data channel. freedv_callback_protorx will be called when a frame
  containing protocol data arrives. freedv_callback_prototx will be called
  when a frame containing protocol information is being generated. Protocol
  information is intended to be used to develop protocols and fancy features
  atop VHF freedv, much like those present in DMR.
   Protocol bits are to be passed in an msb-first char array
   The number of protocol bits are findable with freedv_get_protocol_bits
\*---------------------------------------------------------------------------*/

void freedv_set_callback_protocol(struct freedv *f, freedv_callback_protorx rx, freedv_callback_prototx tx, void *callback_state){
    if (FDV_MODE_ACTIVE( FREEDV_MODE_800XA, f->mode) == false) {
        f->freedv_put_next_proto = rx;
        f->freedv_get_next_proto = tx;
        f->proto_callback_state = callback_state;
    }
}

/*---------------------------------------------------------------------------*\

  FUNCTION....: freedv_set_callback_datarx / freedv_set_callback_datatx
  AUTHOR......: Jeroen Vreeken
  DATE CREATED: 04 March 2016

  Set the callback functions and callback pointer that will be used for the
  data channel. freedv_callback_datarx will be called when a packet has been
  successfully received. freedv_callback_data_tx will be called when 
  transmission of a new packet can begin.
  If the returned size of the datatx callback is zero the data frame is still
  generated, but will contain only a header update.
\*---------------------------------------------------------------------------*/
void freedv_set_callback_data(struct freedv *f, freedv_callback_datarx datarx, freedv_callback_datatx datatx, void *callback_state) {
    if ((FDV_MODE_ACTIVE( FREEDV_MODE_2400A, f->mode)) || (FDV_MODE_ACTIVE( FREEDV_MODE_2400B, f->mode)) || (FDV_MODE_ACTIVE( FREEDV_MODE_800XA, f->mode))){
        if (!f->deframer->fdc)
            f->deframer->fdc = freedv_data_channel_create();
        if (!f->deframer->fdc)
            return;
        
        freedv_data_set_cb_rx(f->deframer->fdc, datarx, callback_state);
        freedv_data_set_cb_tx(f->deframer->fdc, datatx, callback_state);
    }
}

/*---------------------------------------------------------------------------*\

  FUNCTION....: freedv_set_data_header
  AUTHOR......: Jeroen Vreeken
  DATE CREATED: 04 March 2016

  Set the data header for the data channel.
  Header compression will be used whenever packets from this header are sent.
  The header will also be used for fill packets when a data frame is requested
  without a packet available.
\*---------------------------------------------------------------------------*/
void freedv_set_data_header(struct freedv *f, unsigned char *header)
{
    if ((FDV_MODE_ACTIVE( FREEDV_MODE_2400A, f->mode)) || (FDV_MODE_ACTIVE( FREEDV_MODE_2400B, f->mode)) || (FDV_MODE_ACTIVE( FREEDV_MODE_800XA, f->mode))){
        if (!f->deframer->fdc)
            f->deframer->fdc = freedv_data_channel_create();
        if (!f->deframer->fdc)
            return;
        
        freedv_data_set_header(f->deframer->fdc, header);
    }
}

/*---------------------------------------------------------------------------*\

  FUNCTION....: freedv_get_modem_stats
  AUTHOR......: Jim Ahlstrom
  DATE CREATED: 28 July 2015

  Return data from the modem.  The arguments are pointers to the data items.  The
  pointers can be NULL if the data item is not wanted.

\*---------------------------------------------------------------------------*/

void freedv_get_modem_stats(struct freedv *f, int *sync, float *snr_est)
{
    if (FDV_MODE_ACTIVE( FREEDV_MODE_1600, f->mode))
        fdmdv_get_demod_stats(f->fdmdv, &f->stats);
    if ((FDV_MODE_ACTIVE( FREEDV_MODE_700, f->mode)) || (FDV_MODE_ACTIVE( FREEDV_MODE_700B, f->mode))  || (FDV_MODE_ACTIVE( FREEDV_MODE_700C, f->mode)))
        cohpsk_get_demod_stats(f->cohpsk, &f->stats);
    if (FDV_MODE_ACTIVE( FREEDV_MODE_700D, f->mode) || FDV_MODE_ACTIVE( FREEDV_MODE_2020, f->mode)) {
        ofdm_get_demod_stats(f->ofdm, &f->stats);
    }
    if (FDV_MODE_ACTIVE( FREEDV_MODE_2400B, f->mode)) {
        fmfsk_get_demod_stats(f->fmfsk, &f->stats);
    }
    if (sync) *sync = f->stats.sync;
    if (snr_est) *snr_est = f->stats.snr_est;
}

/*---------------------------------------------------------------------------*\

  FUNCTIONS...: freedv_set_*
  AUTHOR......: Jim Ahlstrom
  DATE CREATED: 28 July 2015

  Set some parameters used by FreeDV.  It is possible to write a macro using ## for
  this, but I wasn't sure it would be 100% portable.

\*---------------------------------------------------------------------------*/

// Set integers
void freedv_set_test_frames               (struct freedv *f, int val) {f->test_frames = val;}
void freedv_set_test_frames_diversity	  (struct freedv *f, int val) {f->test_frames_diversity = val;}
void freedv_set_squelch_en                (struct freedv *f, int val) {f->squelch_en = val;}
void freedv_set_total_bit_errors          (struct freedv *f, int val) {f->total_bit_errors = val;}
void freedv_set_total_bits                (struct freedv *f, int val) {f->total_bits = val;}
void freedv_set_total_bit_errors_coded    (struct freedv *f, int val) {f->total_bit_errors_coded = val;}
void freedv_set_total_bits_coded          (struct freedv *f, int val) {f->total_bits_coded = val;}
void freedv_set_clip                      (struct freedv *f, int val) {f->clip = val;}
void freedv_set_varicode_code_num         (struct freedv *f, int val) {varicode_set_code_num(&f->varicode_dec_states, val);}
void freedv_set_ext_vco                   (struct freedv *f, int val) {f->ext_vco = val;}


/* Band Pass Filter to cleanup OFDM tx waveform, only supported by FreeDV 700D */

void freedv_set_tx_bpf(struct freedv *f, int val) {
    if (FDV_MODE_ACTIVE( FREEDV_MODE_700D, f->mode)) {
        ofdm_set_tx_bpf(f->ofdm, val);
    }
}

/* DPSK option for OFDM modem, useful for high SNR, fast fading */

void freedv_set_dpsk(struct freedv *f, int val) {
    if (FDV_MODE_ACTIVE( FREEDV_MODE_700D, f->mode) || FDV_MODE_ACTIVE( FREEDV_MODE_2020, f->mode)) {
        ofdm_set_dpsk(f->ofdm, val);
    }
}

void freedv_set_phase_est_bandwidth_mode(struct freedv *f, int val) {
    if (FDV_MODE_ACTIVE( FREEDV_MODE_700D, f->mode) || FDV_MODE_ACTIVE( FREEDV_MODE_2020, f->mode)) {
        ofdm_set_phase_est_bandwidth_mode(f->ofdm, val);
    }
}

void freedv_set_eq(struct freedv *f, int val) {
    if (FDV_MODE_ACTIVE( FREEDV_MODE_700C, f->mode) || FDV_MODE_ACTIVE( FREEDV_MODE_700D, f->mode)) {
        codec2_700c_eq(f->codec2, val);
    }
}

void freedv_set_verbose(struct freedv *f, int verbosity) {
    f->verbose = verbosity;
    if (FDV_MODE_ACTIVE( FREEDV_MODE_700D, f->mode)) {
        ofdm_set_verbose(f->ofdm, f->verbose);
    }
}

// Set floats
void freedv_set_snr_squelch_thresh        (struct freedv *f, float val) {f->snr_squelch_thresh = val;}

void freedv_set_callback_error_pattern    (struct freedv *f, freedv_calback_error_pattern cb, void *state)
{
    f->freedv_put_error_pattern = cb;
    f->error_pattern_callback_state = state;
}

void freedv_set_carrier_ampl(struct freedv *f, int c, float ampl) {
    assert(FDV_MODE_ACTIVE( FREEDV_MODE_700C, f->mode));
    cohpsk_set_carrier_ampl(f->cohpsk, c, ampl);
}

/*---------------------------------------------------------------------------*\

  FUNCTIONS...: freedv_set_alt_modem_samp_rate
  AUTHOR......: Brady O'Brien
  DATE CREATED: 25 June 2016

  Attempt to set the alternative sample rate on the modem side of the api. Only
   a few alternative sample rates are supported. Please see below.
   
   2400A - 48000, 96000
   2400B - 48000, 96000
  
  TODO: Implement 2400B rate changing, allow other rate changing.
   

\*---------------------------------------------------------------------------*/

int freedv_set_alt_modem_samp_rate(struct freedv *f, int samp_rate){
	if (FDV_MODE_ACTIVE( FREEDV_MODE_2400A, f->mode)){ 
		if(samp_rate == 24000 || samp_rate == 48000 || samp_rate == 96000){
			fsk_destroy(f->fsk);
			f->fsk = fsk_create_hbr(samp_rate,1200,10,4,1200,1200);
        
			FREE(f->tx_bits);
			/* Note: fsk expects tx/rx bits as an array of uint8_ts, not ints */
			f->tx_bits = (int*)MALLOC(f->fsk->Nbits*sizeof(uint8_t));
        
			f->n_nom_modem_samples = f->fsk->N;
			f->n_max_modem_samples = f->fsk->N + (f->fsk->Ts);
			f->n_nat_modem_samples = f->fsk->N;
			f->nin = fsk_nin(f->fsk);
			f->modem_sample_rate = samp_rate;
			return 0;
		}else
			return -1;
	}else if(FDV_MODE_ACTIVE( FREEDV_MODE_2400B, f->mode)){
		if(samp_rate == 48000 || samp_rate == 96000){
			return -1;
		}else
			return -1;
	}
	return -1;
}


/*---------------------------------------------------------------------------* \

  FUNCTIONS...: freedv_set_sync
  AUTHOR......: David Rowe
  DATE CREATED: May 2018

  Extended control of sync state machines, especially for FreeDV 700D.
  This mode is required to acquire sync up at very low SNRS.  This is
  difficult to implement, for example we may get a false sync, or the
  state machine may fall out of sync by mistake during a long fade.

  So with this API call we allow some operator assistance.

  Ensure this is called inthe same thread as freedv_rx().

\*---------------------------------------------------------------------------*/

void freedv_set_sync(struct freedv *freedv, int sync_cmd) {
    assert (freedv != NULL);

    if (FDV_MODE_ACTIVE( FREEDV_MODE_700D, freedv->mode) || FDV_MODE_ACTIVE( FREEDV_MODE_2020, freedv->mode)) {
        ofdm_set_sync(freedv->ofdm, sync_cmd);        
    }
}

struct FSK * freedv_get_fsk(struct freedv *f){
	return f->fsk;
}

/*---------------------------------------------------------------------------*\

  FUNCTIONS...: freedv_get_*
  AUTHOR......: Jim Ahlstrom
  DATE CREATED: 28 July 2015

  Get some parameters from FreeDV.  It is possible to write a macro using ## for
  this, but I wasn't sure it would be 100% portable.

\*---------------------------------------------------------------------------*/

// Get integers
int freedv_get_protocol_bits              (struct freedv *f) {return f->n_protocol_bits;}
int freedv_get_mode                       (struct freedv *f) {return f->mode;}
int freedv_get_test_frames                (struct freedv *f) {return f->test_frames;}
int freedv_get_n_speech_samples           (struct freedv *f) {return f->n_speech_samples;}
int freedv_get_modem_sample_rate          (struct freedv *f) {return f->modem_sample_rate;}
int freedv_get_modem_symbol_rate          (struct freedv *f) {return f->modem_symbol_rate;}
int freedv_get_n_max_modem_samples        (struct freedv *f) {return f->n_max_modem_samples;}
int freedv_get_n_nom_modem_samples        (struct freedv *f) {return f->n_nom_modem_samples;}
int freedv_get_total_bits                 (struct freedv *f) {return f->total_bits;}
int freedv_get_total_bit_errors           (struct freedv *f) {return f->total_bit_errors;}
int freedv_get_total_bits_coded           (struct freedv *f) {return f->total_bits_coded;}
int freedv_get_total_bit_errors_coded     (struct freedv *f) {return f->total_bit_errors_coded;}
int freedv_get_sync                       (struct freedv *f) {return f->stats.sync;}
int freedv_get_speech_sample_rate         (struct freedv *f) {return f->speech_sample_rate;}

int freedv_get_sync_interleaver(struct freedv *f) {
    if (FDV_MODE_ACTIVE( FREEDV_MODE_700D, f->mode) || FDV_MODE_ACTIVE( FREEDV_MODE_2020, f->mode)) {
        return f->ofdm->sync_state_interleaver == synced;
    }

    return 0;
}

int freedv_get_sz_error_pattern(struct freedv *f) 
{
    if (FDV_MODE_ACTIVE( FREEDV_MODE_700, f->mode) || FDV_MODE_ACTIVE( FREEDV_MODE_700B, f->mode) || FDV_MODE_ACTIVE( FREEDV_MODE_700C, f->mode)) {
        /* if diversity disabled callback sends error pattern for upper and lower carriers */
        return f->sz_error_pattern * (2 - f->test_frames_diversity);
    }
    else {
        return f->sz_error_pattern;
    }
}

// Get floats

struct CODEC2 *freedv_get_codec2	(struct freedv *f){return  f->codec2;}
int freedv_get_n_codec_bits             (struct freedv *f){return f->n_codec_bits;}

// Get modem status

void freedv_get_modem_extended_stats(struct freedv *f, struct MODEM_STATS *stats)
{
    if (FDV_MODE_ACTIVE( FREEDV_MODE_1600, f->mode))
        fdmdv_get_demod_stats(f->fdmdv, stats);

    if (FDV_MODE_ACTIVE( FREEDV_MODE_2400A, f->mode) || FDV_MODE_ACTIVE( FREEDV_MODE_800XA, f->mode)) {
        fsk_get_demod_stats(f->fsk, stats);
        float EbNodB = stats->snr_est;                       /* fsk demod actually estimates Eb/No     */
        stats->snr_est = EbNodB + 10.0*log10f(800.0/3000.0); /* so convert to SNR Rb=800, noise B=3000 */
    }

    if (FDV_MODE_ACTIVE( FREEDV_MODE_2400B, f->mode)) {
        fmfsk_get_demod_stats(f->fmfsk, stats);
    }
    
    if (FDV_MODE_ACTIVE( FREEDV_MODE_700, f->mode) || FDV_MODE_ACTIVE( FREEDV_MODE_700B, f->mode) || FDV_MODE_ACTIVE( FREEDV_MODE_700C, f->mode)) {
        cohpsk_get_demod_stats(f->cohpsk, stats);
    }
    
    if (FDV_MODE_ACTIVE( FREEDV_MODE_700D, f->mode) || FDV_MODE_ACTIVE( FREEDV_MODE_2020, f->mode)) {
        ofdm_get_demod_stats(f->ofdm, stats);
    }
    
}

