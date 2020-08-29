/*---------------------------------------------------------------------------*\

  FILE........: freedv_api.c
  AUTHOR......: David Rowe
  DATE CREATED: August 2014

  Library of API functions that implement FreeDV "modes", useful for
  embedding FreeDV in other programs.  Please see:

  1. README_freedv.md
  2. Notes function use in freedv_api.c
  3. The sample freedv_tx.c and freedv_rx.c programs

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

#include "fsk.h"
#include "fmfsk.h"
#include "codec2.h"
#include "codec2_fdmdv.h"
#include "fdmdv_internal.h"
#include "varicode.h"
#include "freedv_api.h"
#include "freedv_api_internal.h"
#include "freedv_vhf_framing.h"
#include "comp_prim.h"

#include "codec2_ofdm.h"
#include "ofdm_internal.h"
#include "mpdecode_core.h"
#include "gp_interleaver.h"
#include "interldpc.h"

#include "debug_alloc.h"

#define VERSION     14    /* The API version number.  The first version
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
 * Version 13   November 2019
 *              Removed 700 and 700B modes
 * Version 14   May 2020
 *              Number of returned speech samples can vary, use freedv_get_n_max_speech_samples() to allocate
 *              buffers.
 */

char *ofdm_statemode[] = {"search","trial","synced"};

char *rx_sync_flags_to_text[] = {
    "----",
    "---T",
    "--S-",
    "--ST",
    "-B--",
    "-B-T",
    "-BS-",
    "-BST",
    "E---",
    "E--T",
    "E-S-",
    "E-ST",
    "EB--",
    "EB-T",
    "EBS-",
    "EBST"};
    
/*---------------------------------------------------------------------------* \

  FUNCTION....: freedv_open
  AUTHOR......: David Rowe
  DATE CREATED: 3 August 2014

  Call this first to initialise.  Returns NULL if initialisation
  fails. If a malloc() or calloc() fails in general asserts() will
  fire.

\*---------------------------------------------------------------------------*/

struct freedv *freedv_open(int mode) {
    return freedv_open_advanced(mode, NULL);
}

// Ignoring adv parameter
struct freedv *freedv_open_advanced(int mode, struct freedv_advanced *adv) {
    struct freedv *f;
    
    if (false == (FDV_MODE_ACTIVE( FREEDV_MODE_1600,mode) || FDV_MODE_ACTIVE( FREEDV_MODE_2400A,mode) || 
		FDV_MODE_ACTIVE( FREEDV_MODE_2400B,mode) || FDV_MODE_ACTIVE( FREEDV_MODE_800XA,mode) || 
		FDV_MODE_ACTIVE( FREEDV_MODE_700C,mode) || FDV_MODE_ACTIVE( FREEDV_MODE_700D,mode)  ||
                  FDV_MODE_ACTIVE( FREEDV_MODE_2020,mode)) ) return NULL;

    /* set everything to zero just in case */
    f = (struct freedv*)CALLOC(1, sizeof(struct freedv));
    if (f == NULL) return NULL;

    f->mode = mode;

    if (FDV_MODE_ACTIVE( FREEDV_MODE_1600, mode)) freedv_1600_open(f);
    if (FDV_MODE_ACTIVE( FREEDV_MODE_700C, mode)) freedv_700c_open(f);
    if (FDV_MODE_ACTIVE( FREEDV_MODE_700D, mode)) freedv_700d_open(f);
#ifdef __LPCNET__
    if (FDV_MODE_ACTIVE( FREEDV_MODE_2020, mode)) freedv_2020_open(f);
#endif    
    if (FDV_MODE_ACTIVE( FREEDV_MODE_2400A, mode)) freedv_2400a_open(f);
    if (FDV_MODE_ACTIVE( FREEDV_MODE_2400B, mode) ) freedv_2400b_open(f);
    if (FDV_MODE_ACTIVE( FREEDV_MODE_800XA, mode)) freedv_800xa_open(f);

    varicode_decode_init(&f->varicode_dec_states, 1);

    return f;
}

/*---------------------------------------------------------------------------*\

  FUNCTION....: freedv_close
  AUTHOR......: David Rowe
  DATE CREATED: 3 August 2014

  Call to shut down a freedv instance and free memory.

\*---------------------------------------------------------------------------*/

void freedv_close(struct freedv *freedv) {
    assert(freedv != NULL);

    FREE(freedv->tx_payload_bits);
    FREE(freedv->rx_payload_bits);
    if (freedv->codec2) codec2_destroy(freedv->codec2);

    if (FDV_MODE_ACTIVE(FREEDV_MODE_1600, freedv->mode)) {
        FREE(freedv->fdmdv_bits);
        FREE(freedv->fdmdv_tx_bits);
        FREE(freedv->fdmdv_rx_bits);
        fdmdv_destroy(freedv->fdmdv);
    }

    if (FDV_MODE_ACTIVE( FREEDV_MODE_700C, freedv->mode)) {
        cohpsk_destroy(freedv->cohpsk);
        quisk_filt_destroy(freedv->ptFilter8000to7500);
        FREE(freedv->ptFilter8000to7500);
        quisk_filt_destroy(freedv->ptFilter7500to8000);
        FREE(freedv->ptFilter7500to8000);
    }

    if (FDV_MODE_ACTIVE( FREEDV_MODE_700D, freedv->mode)) {
        FREE(freedv->codeword_symbols);
        FREE(freedv->codeword_amps);
        FREE(freedv->ldpc);
        ofdm_destroy(freedv->ofdm);
    }

    if (FDV_MODE_ACTIVE( FREEDV_MODE_2020, freedv->mode)) {
        FREE(freedv->codeword_symbols);
        FREE(freedv->codeword_amps);
        FREE(freedv->ldpc);
        FREE(freedv->passthrough_2020);
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
    
    FREE(freedv);
}


/* helper function, unpacked bits are much easier to work with inside the modem */

static void codec2_encode_upacked(struct freedv *f, uint8_t unpacked_bits[], short speech_in[]) {
    int n_packed = (f->bits_per_codec_frame + 7) / 8;
#ifdef NO_C99
    uint8_t *packed_codec_bits = alloca(n_packed);
#else
    uint8_t packed_codec_bits[n_packed];
#endif
    
    codec2_encode(f->codec2, packed_codec_bits, speech_in);
    
    /* unpack bits, MSB first */
    int bit = 7, byte = 0;
    for(int i=0; i<f->bits_per_codec_frame; i++) {
        unpacked_bits[i] = (packed_codec_bits[byte] >> bit) & 0x1;
        bit--;
        if (bit < 0) {
            bit = 7;
            byte++;
        }
    }
}


/*---------------------------------------------------------------------------*\

  FUNCTION....: freedv_tx
  AUTHOR......: David Rowe
  DATE CREATED: 3 August 2014

  Takes a frame of input speech samples, encodes and modulates them to
  produce a frame of modem samples that can be sent to the
  transmitter.  See freedv_tx.c for an example.

  speech_in[] is sampled at freedv_get_speech_sample_rate() Hz, and the
  user must supply a block of exactly
  freedv_get_n_speech_samples(). The speech_in[] level should be such
  that the peak speech level is between +/- 16384 and +/- 32767.

  The FDM modem signal mod_out[] is sampled at
  freedv_get_modem_sample_rate() and is always exactly
  freedv_get_n_nom_modem_samples() long.  mod_out[] will be scaled
  such that the peak level is just less than +/-32767.

  The FreeDV 1600/700C/700D/2020 waveforms have a crest factor of
  around 10dB, similar to SSB.  These modes are usually operated at a
  "backoff" of 8dB.  Adjust the power amplifier drive so that the
  average power is 8dB less than the peak power of the PA.  For
  example, on a radio rated at 100W PEP for SSB, the average FreeDV
  power is typically 20W.

  Caution - some PAs cannot handle a high continuous power.  A
  conservative level is 20W average for a 100W PEP rated PA.

  The FreeDV 2400A/800XA modes are constant amplitude, designed for
  Class C PAs.  They have a crest factor of 3dB. If using a SSB PA,
  adjust the drive so you average power is within the limits of your PA
  (e.g. 20W average for a 100W PA).

\*---------------------------------------------------------------------------*/

/* real-valued short output */

void freedv_tx(struct freedv *f, short mod_out[], short speech_in[]) {
    assert(f != NULL);
#ifdef NO_C99
    COMP *tx_fdm = alloca(f->n_nom_modem_samples*sizeof(COMP));
#else
    COMP tx_fdm[f->n_nom_modem_samples];
#endif
    int  i;
    assert((FDV_MODE_ACTIVE( FREEDV_MODE_1600, f->mode))  || (FDV_MODE_ACTIVE( FREEDV_MODE_700C, f->mode))  ||
           (FDV_MODE_ACTIVE( FREEDV_MODE_700D, f->mode))  || 
           (FDV_MODE_ACTIVE( FREEDV_MODE_2400A, f->mode)) || (FDV_MODE_ACTIVE( FREEDV_MODE_2400B, f->mode)) || 
           (FDV_MODE_ACTIVE( FREEDV_MODE_800XA, f->mode)) || (FDV_MODE_ACTIVE( FREEDV_MODE_2020, f->mode)));
    
    /* FSK and MEFSK/FMFSK modems work only on real samples. It's simpler to just 
     * stick them in the real sample tx/rx functions than to add a comp->real converter
     * to comptx */
     
    if ((FDV_MODE_ACTIVE( FREEDV_MODE_2400A, f->mode)) || (FDV_MODE_ACTIVE( FREEDV_MODE_2400B, f->mode)) || (FDV_MODE_ACTIVE( FREEDV_MODE_800XA, f->mode))){
        /* 800XA has two codec frames per modem frame */
        if(FDV_MODE_ACTIVE( FREEDV_MODE_800XA, f->mode)){
            codec2_encode(f->codec2, &f->tx_payload_bits[0], &speech_in[  0]);
            codec2_encode(f->codec2, &f->tx_payload_bits[4], &speech_in[320]);
        }else{
            codec2_encode(f->codec2, f->tx_payload_bits, speech_in);
        }
        freedv_tx_fsk_voice(f, mod_out);
    } else {
        freedv_comptx(f, tx_fdm, speech_in);
        for(i=0; i<f->n_nom_modem_samples; i++)
            mod_out[i] = tx_fdm[i].real;
    }
}


/* complex float output samples version */

void freedv_comptx(struct freedv *f, COMP mod_out[], short speech_in[]) {
    assert(f != NULL);

    assert((FDV_MODE_ACTIVE( FREEDV_MODE_1600, f->mode)) || (FDV_MODE_ACTIVE( FREEDV_MODE_700C, f->mode)) || 
           (FDV_MODE_ACTIVE( FREEDV_MODE_2400A, f->mode)) || (FDV_MODE_ACTIVE( FREEDV_MODE_2400B, f->mode)) ||
           (FDV_MODE_ACTIVE( FREEDV_MODE_700D, f->mode))  || (FDV_MODE_ACTIVE( FREEDV_MODE_2020, f->mode)));

    if (FDV_MODE_ACTIVE( FREEDV_MODE_1600, f->mode)) {
        codec2_encode_upacked(f, f->tx_payload_bits, speech_in);
        freedv_comptx_fdmdv_1600(f, mod_out);
    }
    
    /* all these modes need to pack a bunch of codec frames into one modem frame ... */
    
    if (FDV_MODE_ACTIVE( FREEDV_MODE_700C, f->mode)) {
        for (int j=0; j<f->n_codec_frames; j++) {
            codec2_encode_upacked(f, f->tx_payload_bits+j*f->bits_per_codec_frame, speech_in);
            speech_in += codec2_samples_per_frame(f->codec2);
        }
        freedv_comptx_700c(f, mod_out);
    }

    if (FDV_MODE_ACTIVE( FREEDV_MODE_700D, f->mode)) {

        /* buffer up bits until we get enough encoded bits for interleaver */
        
        for (int j=0; j<f->n_codec_frames; j++) {
            int offset = (f->modem_frame_count_tx*f->n_codec_frames + j)*f->bits_per_codec_frame;
            codec2_encode_upacked(f, f->tx_payload_bits + offset, speech_in);
            speech_in += codec2_samples_per_frame(f->codec2);
        }

        freedv_comptx_700d(f, mod_out);
    }
    
#ifdef __LPCNET__
    if (FDV_MODE_ACTIVE( FREEDV_MODE_2020, f->mode)) {
       
        /* buffer up bits until we get enough encoded bits for interleaver */
        
        for (int j=0; j<f->n_codec_frames; j++) {
            int offset = (f->modem_frame_count_tx*f->n_codec_frames + j)*f->bits_per_codec_frame;
            lpcnet_enc(f->lpcnet, speech_in, (char*)f->tx_payload_bits + offset);
            speech_in += lpcnet_samples_per_frame(f->lpcnet);
        }

        freedv_comptx_2020(f, mod_out);
    }
#endif
    
    /* 2400 A and B are handled by the real-mode TX */
    if(FDV_MODE_ACTIVE( FREEDV_MODE_2400A, f->mode) || FDV_MODE_ACTIVE( FREEDV_MODE_2400B, f->mode)){
        codec2_encode(f->codec2, f->tx_payload_bits, speech_in);
        freedv_comptx_fsk_voice(f, mod_out);
    }
}

/* a way to send raw frames of bytes, or speech data that was compressed externally, real short output */
void freedv_rawdatatx(struct freedv *f, short mod_out[], unsigned char *packed_payload_bits) {
    assert(f != NULL);
#ifdef NO_C99
    COMP *tx_fdm = alloca(f->n_nom_modem_samples*sizeof(COMP));
#else
    COMP tx_fdm[f->n_nom_modem_samples];
#endif

    /* FSK modes used packed bits */
    if(FDV_MODE_ACTIVE( FREEDV_MODE_2400A, f->mode) || FDV_MODE_ACTIVE( FREEDV_MODE_2400B, f->mode) ||
       FDV_MODE_ACTIVE( FREEDV_MODE_800XA, f->mode) ) {
	freedv_codec_frames_from_rawdata(f, f->tx_payload_bits,  packed_payload_bits);
        freedv_tx_fsk_voice(f, mod_out);
        return; /* output is already real */
    }

    /* unpack bits, MSB first */
    int bit = 7, byte = 0;
    for(int i=0; i<f->bits_per_modem_frame; i++) {
        f->tx_payload_bits[i] = (packed_payload_bits[byte] >> bit) & 0x1;
        bit--;
        if (bit < 0) {
            bit = 7;
            byte++;
        }
    }

    if (FDV_MODE_ACTIVE( FREEDV_MODE_1600, f->mode)) freedv_comptx_fdmdv_1600(f, tx_fdm);
    if (FDV_MODE_ACTIVE( FREEDV_MODE_700C, f->mode)) freedv_comptx_700c(f, tx_fdm);
    if (FDV_MODE_ACTIVE( FREEDV_MODE_700D, f->mode)) freedv_comptx_700d(f, tx_fdm);

    /* convert complex to real */
    for(int i=0; i<f->n_nom_modem_samples; i++)
        mod_out[i] = tx_fdm[i].real;
}


/* VHF packet data tx function */
void freedv_datatx (struct freedv *f, short mod_out[]) {
    assert(f != NULL);
    if (FDV_MODE_ACTIVE( FREEDV_MODE_2400A, f->mode) || FDV_MODE_ACTIVE( FREEDV_MODE_2400B, f->mode) || FDV_MODE_ACTIVE( FREEDV_MODE_800XA, f->mode)) {
        freedv_tx_fsk_data(f, mod_out);
    }
}


/* VHF packet data: returns how many tx frames are queued up but not sent yet */
int  freedv_data_ntxframes (struct freedv *f) {
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
    if (FDV_MODE_ACTIVE( FREEDV_MODE_700C, f->mode))
        // For mode 700C, the input rate is 8000 sps, but the modem rate is 7500 sps
        // For mode 700C, we request a larger number of Rx samples that will be decimated to f->nin samples
        return (16 * f->nin + f->ptFilter8000to7500->decim_index) / 15;
    else
        return f->nin;
}

int freedv_codec_frames_from_rawdata(struct freedv *f, unsigned char *codec_frames, unsigned char *rawdata)
{
    int cbit = 7;
    int cbyte = 0;
    int rbit = 7;
    int rbyte = 0;
    int modem_bits = freedv_get_bits_per_modem_frame(f);
    int codec_bits = freedv_get_bits_per_codec_frame(f);
    int nr_cbits = 0;
    int i;
    
    codec_frames[0] = 0;
    for (i = 0; i < modem_bits; i++) {
        codec_frames[cbyte] |= ((rawdata[rbyte] >> rbit) & 1) << cbit;

        rbit--;
        if (rbit < 0) {
            rbit = 7;
            rbyte++;
        }
	
        cbit--;
        if (cbit < 0) {
            cbit = 7;
	    cbyte++;
	    codec_frames[cbyte] = 0;
	}
	nr_cbits++;
	if (nr_cbits == codec_bits) {
            if (cbit) {
                cbyte++;
                codec_frames[cbyte] = 0;
            }
            cbit = 7;
	    nr_cbits = 0;
	}
    }
    return f->n_codec_frames;
}

int freedv_rawdata_from_codec_frames(struct freedv *f, unsigned char *rawdata, unsigned char *codec_frames)
{
    int cbit = 7;
    int cbyte = 0;
    int rbit = 7;
    int rbyte = 0;
    int modem_bits = freedv_get_bits_per_modem_frame(f);
    int codec_bits = freedv_get_bits_per_codec_frame(f);
    int nr_cbits = 0;
    int i;
    
    rawdata[rbyte] = 0;
    for (i = 0; i < modem_bits; i++) {
        rawdata[rbyte] |= ((codec_frames[cbyte] >> cbit) & 1) << rbit;

        rbit--;
        if (rbit < 0) {
            rbit = 7;
            rbyte++;
	    rawdata[rbyte] = 0;
        }
	
        cbit--;
        if (cbit < 0) {
            cbit = 7;
	    cbyte++;
	}
	
	nr_cbits++;
	if (nr_cbits == codec_bits) {
            if (cbit)
                cbyte++;
            cbit = 7;
	    nr_cbits = 0;
	}
    }
    return f->n_codec_frames;
}


/*---------------------------------------------------------------------------*\

  FUNCTION....: freedv_rx
  AUTHOR......: David Rowe
  DATE CREATED: 3 August 2014

  Takes samples from the radio receiver, demodulates and FEC decodes
  them, producing a frame of decoded speech samples.  See freedv_rx.c
  for an example.

  demod_in[] is a block of received samples sampled at
  freedv_get_modem_sample_rate().  To account for difference in the
  transmit and receive sample clock frequencies, the number of
  demod_in[] samples is time varying. You MUST call freedv_nin()
  BEFORE EACH call to freedv_rx() and pass exactly that many samples
  to this function:

  short demod_in[freedv_get_n_max_modem_samples(f)];
  short speech_out[freedv_get_n_max_speech_samples(f)];

  nin = freedv_nin(f);
  while(fread(demod_in, sizeof(short), nin, fin) == nin) {
      nout = freedv_rx(f, speech_out, demod_in);
      fwrite(speech_out, sizeof(short), nout, fout);
      nin = freedv_nin(f);
  }

  To help set your buffer sizes, The maximum value of freedv_nin() is
  freedv_get_n_max_modem_samples().

  freedv_rx() returns the number of output speech samples available in
  speech_out[], which is sampled at freedv_get_speech_sample_rate(f).
  You should ALWAYS check the return value of freedv_rx(), and read
  EXACTLY that number of speech samples from speech_out[].  

  Not every call to freedv_rx will return speech samples; in some
  modes several modem frames are processed before speech samples are
  returned.  When squelch is active, zero samples may be returned.

  1600 and 700D mode: When out of sync, the number of output speech
  samples returned will be freedv_nin(). When in sync to a valid
  FreeDV 1600 signal, the number of output speech samples will
  alternate between freedv_get_n_speech_samples() and 0.

  The peak level of demod_in[] is not critical, as the demod works
  well over a wide range of amplitude scaling.  However avoid clipping
  (overload, or samples pinned to +/- 32767).  speech_out[] will peak
  at just less than +/-32767.

  When squelch is disabled, this function echoes the demod_in[]
  samples to speech_out[].  This allows the user to listen to the
  channel, which is useful for tuning FreeDV signals or reception of
  non-FreeDV signals.

\*---------------------------------------------------------------------------*/

int freedv_rx(struct freedv *f, short speech_out[], short demod_in[]) {
    assert(f != NULL);
    int i;
    int nin = freedv_nin(f);
    f->nin_prev = nin;

    assert(nin <= f->n_max_modem_samples);
       
    /* FSK RX happens in real floats, so convert to those and call their demod here */
    if( (FDV_MODE_ACTIVE( FREEDV_MODE_2400A, f->mode)) || (FDV_MODE_ACTIVE( FREEDV_MODE_2400B, f->mode)) || (FDV_MODE_ACTIVE( FREEDV_MODE_800XA, f->mode)) ){
#ifdef NO_C99
        float *rx_float = alloca(f->n_max_modem_samples*sizeof(float));
#else
        float rx_float[f->n_max_modem_samples];
#endif
        for(i=0; i<nin; i++) {
            rx_float[i] = ((float)demod_in[i]);
        }
        return freedv_floatrx(f,speech_out,rx_float);
    }
    
    if ( (FDV_MODE_ACTIVE( FREEDV_MODE_1600, f->mode)) || (FDV_MODE_ACTIVE( FREEDV_MODE_700C, f->mode))
         || (FDV_MODE_ACTIVE( FREEDV_MODE_2020, f->mode))) {

        float gain = 1.0f;
        
        assert(nin <= f->n_max_modem_samples);
#ifdef NO_C99
        COMP *rx_fdm = alloca(f->n_max_modem_samples*sizeof(COMP));
#else
        COMP rx_fdm[f->n_max_modem_samples];
#endif
        
        for(i=0; i<nin; i++) {
            rx_fdm[i].real = gain*(float)demod_in[i];
            rx_fdm[i].imag = 0.0f;
        }
        return freedv_comprx(f, speech_out, rx_fdm);
    }

    /* special low memory version for 700D, to help with stm32 port */
    if (FDV_MODE_ACTIVE( FREEDV_MODE_700D, f->mode)) {
        float gain = 2.0f; /* keep levels the same as Octave simulations and C unit tests for real signals */
        return freedv_shortrx(f, speech_out, demod_in, gain);
    }
    
    assert(1); /* should never get here */
    return 0;
}

/* complex sample input version from the radio */

int freedv_comprx(struct freedv *f, short speech_out[], COMP demod_in[]) {
    assert(f != NULL);    
    assert(f->nin <= f->n_max_modem_samples);
    int rx_status = 0;
    f->nin_prev = freedv_nin(f);
   
    if (FDV_MODE_ACTIVE( FREEDV_MODE_1600, f->mode)) {
        rx_status = freedv_comprx_fdmdv_1600(f, demod_in);
    }
    if (FDV_MODE_ACTIVE( FREEDV_MODE_700C, f->mode)) {
        rx_status = freedv_comprx_700c(f, demod_in);
    }

    if( (FDV_MODE_ACTIVE( FREEDV_MODE_2400A, f->mode)) || (FDV_MODE_ACTIVE( FREEDV_MODE_2400B, f->mode)) || (FDV_MODE_ACTIVE( FREEDV_MODE_800XA, f->mode))){
        rx_status = freedv_comprx_fsk(f, demod_in);
    }

    if (FDV_MODE_ACTIVE( FREEDV_MODE_700D, f->mode)) {
        rx_status = freedv_comp_short_rx_700d(f, (void*)demod_in, 0, 2.0f); // was 1.0 ??
    }
    
    if (FDV_MODE_ACTIVE( FREEDV_MODE_2020, f->mode)) {
#ifdef __LPCNET__
        rx_status = freedv_comprx_2020(f, demod_in);
#endif
    }

#ifdef NO_C99
    short *demod_in_short = alloca(f->nin_prev*sizeof(short));
#else
    short demod_in_short[f->nin_prev];
#endif
    
    for(int i=0; i<f->nin_prev; i++)
        demod_in_short[i] = demod_in[i].real;
    
    return freedv_bits_to_speech(f, speech_out, demod_in_short, rx_status);
}

/* memory efficient real short version - just for 700D on the SM1000 */

int freedv_shortrx(struct freedv *f, short speech_out[], short demod_in[], float gain) {
    assert(f != NULL);
    int rx_status;
    f->nin_prev = f->nin;

    // At this stage short interface only supported for 700D, to help
    // memory requirements on stm32
    assert(f->mode == FREEDV_MODE_700D);
    assert(f->nin <= f->n_max_modem_samples);

    if (FDV_MODE_ACTIVE( FREEDV_MODE_700D, f->mode)) {
        rx_status = freedv_comp_short_rx_700d(f, (void*)demod_in, 1, gain);
    }
    
    return freedv_bits_to_speech(f, speech_out, demod_in, rx_status);
}


/* helper function, unpacked bits are much easier to work with inside the modem */

static void codec2_decode_upacked(struct freedv *f, short speech_out[], uint8_t unpacked_bits[]) {
    int n_packed = (f->bits_per_codec_frame + 7) / 8;
#ifdef NO_C99
    uint8_t *packed_codec_bits = alloca(n_packed*sizeof(uint8_t));
#else
    uint8_t packed_codec_bits[n_packed];
#endif
    memset(packed_codec_bits, 0, n_packed);
    
    /* pack bits, MSB received first */
    int bit = 7, byte = 0;
    for(int i=0; i<f->bits_per_codec_frame; i++) {
        packed_codec_bits[byte] |= (unpacked_bits[i] << bit);
        bit--;
        if (bit < 0) {
            bit = 7;
            byte++;
        }
    }
    codec2_decode(f->codec2, speech_out, packed_codec_bits);    
}


/*---------------------------------------------------------------------------* \

  FUNCTION....: freedv_rx_bits_to_speech
  AUTHOR......: David Rowe
  DATE CREATED: May 2020

  The *_rx functions takes off air samples, demodulate and (for some
  modes) FEC decode, giving us a frame of bits.

  This function captures a lot of tricky logic that has been distilled
  through experience:

  There may not be a frame of bits returned on every call freedv_*rx* call.
  When there are valid bits we need to run the speech decoder.
  We may not have demod sync, so various pass through options may happen
  with the input samples.  
  We may squelch based on SNR.
  Need to handle various codecs, and varying number of codec frames per modem frame
  Squelch audio if test frames are being sent
  Determine how many speech samples to return, which will vary if in sync/out of sync
  Work with real and complex inputs (complex wrapper)
  Attenuate audio on pass through
  Deal with 700D first frame burble, and different sync states from OFDM modes like 700D
  Output no samples if squelched, we assume it's OK for the audio sink to run dry
  A FIFO required on output to smooth sample flow to audio sink
  Don't decode when we are sendingtest frames

\*---------------------------------------------------------------------------*/

int freedv_bits_to_speech(struct freedv *f, short speech_out[], short demod_in[], int rx_status) {
    int nout = 0;
    int decode_speech = 0;
   
    if ((rx_status & RX_SYNC) == 0) {

        if (f->squelch_en == 0) {

            /* attenuate audio 12dB bit as channel noise isn't that pleasant */
            float passthrough_gain = 0.25;
            
            /* pass through received samples so we can hear what's going on, e.g. during tuning */

            if (f->mode == FREEDV_MODE_2020) {           
                /* 8kHz modem sample rate but 16 kHz speech sample
                   rate, so we need to resample */
                nout = 2*f->nin_prev;
                assert(nout <= freedv_get_n_max_speech_samples(f));
#ifdef NO_C99
                float *tmp = alloca(nout*sizeof(float));
#else
                float tmp[nout];
#endif
                for(int i=0; i<nout/2; i++)
                    f->passthrough_2020[FDMDV_OS_TAPS_16K+i] = demod_in[i];
                fdmdv_8_to_16(tmp, &f->passthrough_2020[FDMDV_OS_TAPS_16K], nout/2);
                for(int i=0; i<nout; i++)
                    speech_out[i] = passthrough_gain*tmp[i];
            } else {
	        /* Speech and modem rates might be different */
	        int rate_factor = f->modem_sample_rate / f-> speech_sample_rate;
                nout = f->nin_prev / rate_factor;
                for(int i=0; i<nout; i++)
                    speech_out[i] = passthrough_gain*demod_in[i * rate_factor];
           }
        }
    }

    if ((rx_status & RX_SYNC) && (rx_status & RX_BITS) && !f->test_frames) {
       /* following logic is tricky so spell it out clearly, see table
          in: https://github.com/drowe67/codec2/pull/111 */
        
       if (f->squelch_en == 0) {
           decode_speech = 1;
       } else {
           /* squelch is enabled */
           
           /* anti-burble case - don't decode on trial sync unless the
              frame has no bit errors.  This prevents short lived trial
              sync cases generating random bursts of audio */
           if (rx_status & RX_TRIAL_SYNC) {
               if ((rx_status & RX_BIT_ERRORS) == 0)
                   decode_speech = 1;
           }
           else {
               /* sync is solid - decode even through fades as there is still some speech info there */
               if (f->stats.snr_est > f->snr_squelch_thresh)
                   decode_speech = 1;
           }
       }
    }
   
    if (decode_speech) {
        if(FDV_MODE_ACTIVE( FREEDV_MODE_2020, f->mode)) {
#ifdef __LPCNET__
            /* LPCNet decoder */
            
            int bits_per_codec_frame = lpcnet_bits_per_frame(f->lpcnet);
            int data_bits_per_frame = f->ldpc->data_bits_per_frame;
            int frames = data_bits_per_frame/bits_per_codec_frame;            

            nout = f->n_speech_samples;
            for (int i = 0; i < frames; i++) {
                lpcnet_dec(f->lpcnet, (char*) f->rx_payload_bits + (i + frames * f->modem_frame_count_rx) * bits_per_codec_frame, speech_out);
                speech_out += lpcnet_samples_per_frame(f->lpcnet);
            }
            f->modem_frame_count_rx++;
#endif    
        } else {
            /* codec 2 decoder */

            if(FDV_MODE_ACTIVE( FREEDV_MODE_700D, f->mode)) {
                nout = f->n_speech_samples;
                for (int i = 0; i < f->n_codec_frames; i++) {
                    codec2_decode_upacked(f, speech_out, f->rx_payload_bits + (i + f->n_codec_frames * f->modem_frame_count_rx) * f->bits_per_codec_frame);
                    speech_out += codec2_samples_per_frame(f->codec2);
                }
                f->modem_frame_count_rx++;
            } else { 
                /* non-interleaved Codec 2 modes */
                
                nout = f->n_speech_samples;
                if ( (FDV_MODE_ACTIVE( FREEDV_MODE_2400A, f->mode)) || (FDV_MODE_ACTIVE( FREEDV_MODE_2400B, f->mode)))
                    codec2_decode(f->codec2, speech_out, f->rx_payload_bits);
                else if (FDV_MODE_ACTIVE( FREEDV_MODE_800XA, f->mode)) {
                    codec2_decode(f->codec2, &speech_out[  0], &f->rx_payload_bits[0]);
                    codec2_decode(f->codec2, &speech_out[320], &f->rx_payload_bits[4]);
               } else {
                    for (int i = 0; i <f->n_codec_frames; i++) {
                        codec2_decode_upacked(f, speech_out, f->rx_payload_bits + i*f->bits_per_codec_frame);
                        speech_out += codec2_samples_per_frame(f->codec2);
                    }
                }
            }
        }
    }

    if (f->verbose == 2) {
        fprintf(stderr, "    sqen: %d nout: %d decsp: %d\n", f->squelch_en, nout, decode_speech);
    }
    
    f->rx_status= rx_status;
    assert(nout <= freedv_get_n_max_speech_samples(f));    
    return nout;
}


/* a way to receive raw frames of bytes, or speech data that will be decompressed externally */
int freedv_rawdatarx(struct freedv *f, unsigned char *packed_payload_bits, short demod_in[])
{
    assert(f != NULL);
    int i;
    int nin = freedv_nin(f);
    int ret = 0;
    int rx_status = 0;
    
    assert(nin <= f->n_max_modem_samples);
    f->nin_prev = nin;
    
#ifdef NO_C99
    COMP *rx_fdm = alloca(f->n_max_modem_samples*sizeof(COMP));
#else
    COMP rx_fdm[f->n_max_modem_samples];
#endif
    
    for(i=0; i<nin; i++) {
        rx_fdm[i].real = (float)demod_in[i];
        rx_fdm[i].imag = 0.0;
    }

    /* FSK modes used packed bits internally */
    if (FDV_MODE_ACTIVE( FREEDV_MODE_2400A, f->mode) || FDV_MODE_ACTIVE( FREEDV_MODE_2400B, f->mode) || FDV_MODE_ACTIVE( FREEDV_MODE_800XA, f->mode)){
        rx_status = freedv_comprx_fsk(f, rx_fdm);
        f->rx_status = rx_status;
        if (rx_status & RX_BITS) {
	    ret = (freedv_get_bits_per_modem_frame(f) + 7) / 8;
	    freedv_rawdata_from_codec_frames(f, packed_payload_bits, f->rx_payload_bits);
        }
        return ret;
    }

    if (FDV_MODE_ACTIVE( FREEDV_MODE_1600, f->mode)) rx_status = freedv_comprx_fdmdv_1600(f, rx_fdm);
    if (FDV_MODE_ACTIVE( FREEDV_MODE_700C, f->mode)) rx_status = freedv_comprx_700c(f, rx_fdm);
    if (FDV_MODE_ACTIVE( FREEDV_MODE_700D, f->mode)) rx_status = freedv_comp_short_rx_700d(f, (void*)demod_in, 1, 2.0f);

    if (rx_status & RX_BITS) {
	ret = (f->bits_per_modem_frame+7)/8;
        memset(packed_payload_bits, 0, ret);
        /* pack bits */
        int bit = 7, byte = 0;
        for(int i=0; i<f->bits_per_modem_frame; i++) {
            packed_payload_bits[byte] |= (f->rx_payload_bits[i] << bit);
            bit--;
            if (bit < 0) {
                bit = 7;
                byte++;
            }
        }
    }

    /* might want to check this for errors, e.g. if reliable data is important */
    f->rx_status= rx_status;
    
    return ret;
}


/*---------------------------------------------------------------------------* \

  FUNCTION....: freedv_get_version
  AUTHOR......: Jim Ahlstrom
  DATE CREATED: 28 July 2015

  Return the version of the FreeDV API.  This is meant to help API
  users determine when incompatible changes have occurred.

\*---------------------------------------------------------------------------*/

int freedv_get_version(void)
{
    return VERSION;
}

/*---------------------------------------------------------------------------* \

  FUNCTION....: freedv_get_hash
  AUTHOR......: David Rowe
  DATE CREATED: July 2020

  Return the a string with the Git hash of the repo used to build this code.

\*---------------------------------------------------------------------------*/

static char git_hash[] = GIT_HASH;
char *freedv_get_hash(void)
{
    return git_hash;
}

/*---------------------------------------------------------------------------*\

  FUNCTION....: freedv_set_callback_txt
  AUTHOR......: Jim Ahlstrom
  DATE CREATED: 28 July 2015

  Set the callback functions and the callback state pointer that will
  be used for the aux txt channel.  The freedv_callback_rx is a
  function pointer that will be called to return received characters.
  The freedv_callback_tx is a function pointer that will be called to
  send transmitted characters.  The callback state is a user-defined
  void pointer that will be passed to the callback functions.  Any or
  all can be NULL, and the default is all NULL.

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

  VHF packet data function.

  Set the callback functions and callback pointer that will be used
  for the protocol data channel. freedv_callback_protorx will be
  called when a frame containing protocol data
  arrives. freedv_callback_prototx will be called when a frame
  containing protocol information is being generated. Protocol
  information is intended to be used to develop protocols and fancy
  features atop VHF freedv, much like those present in DMR.  Protocol
  bits are to be passed in an msb-first char array The number of
  protocol bits are findable with freedv_get_protocol_bits

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

  VHF packet data function.

  Set the callback functions and callback pointer that will be used
  for the data channel. freedv_callback_datarx will be called when a
  packet has been successfully received. freedv_callback_data_tx will
  be called when transmission of a new packet can begin.  If the
  returned size of the datatx callback is zero the data frame is still
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

  VHF packet data function.

  Set the data header for the data channel.  Header compression will
  be used whenever packets from this header are sent.  The header will
  also be used for fill packets when a data frame is requested without
  a packet available.

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

  Return data from the modem.  The arguments are pointers to the data
  items.  The pointers can be NULL if the data item is not wanted.

\*---------------------------------------------------------------------------*/

void freedv_get_modem_stats(struct freedv *f, int *sync, float *snr_est)
{
    if (FDV_MODE_ACTIVE( FREEDV_MODE_1600, f->mode))
        fdmdv_get_demod_stats(f->fdmdv, &f->stats);
    if (FDV_MODE_ACTIVE( FREEDV_MODE_700C, f->mode))
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

  Set some parameters used by FreeDV.  It is possible to write a macro
  using ## for this, but I wasn't sure it would be 100% portable.

\*---------------------------------------------------------------------------*/

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
void freedv_set_snr_squelch_thresh        (struct freedv *f, float val) {f->snr_squelch_thresh = val;}

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
    if (FDV_MODE_ACTIVE( FREEDV_MODE_700C, f->mode)) {
        cohpsk_set_verbose(f->cohpsk, f->verbose);
    }
    if (FDV_MODE_ACTIVE( FREEDV_MODE_700D, f->mode)) {
        ofdm_set_verbose(f->ofdm, f->verbose);
    }
}

void freedv_set_callback_error_pattern(struct freedv *f, freedv_calback_error_pattern cb, void *state)
{
    f->freedv_put_error_pattern = cb;
    f->error_pattern_callback_state = state;
}

void freedv_set_carrier_ampl(struct freedv *f, int c, float ampl) {
    assert(FDV_MODE_ACTIVE( FREEDV_MODE_700C, f->mode));
    cohpsk_set_carrier_ampl(f->cohpsk, c, ampl);
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

  Ensure this is called in the same thread as freedv_rx().

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

  Get some parameters from FreeDV.

\*---------------------------------------------------------------------------*/

int freedv_get_protocol_bits              (struct freedv *f) {return f->n_protocol_bits;}
int freedv_get_mode                       (struct freedv *f) {return f->mode;}
int freedv_get_test_frames                (struct freedv *f) {return f->test_frames;}
int freedv_get_speech_sample_rate         (struct freedv *f) {return f-> speech_sample_rate;}
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
struct CODEC2 *freedv_get_codec2	  (struct freedv *f){return  f->codec2;}
int freedv_get_bits_per_codec_frame       (struct freedv *f){return f->bits_per_codec_frame;}
int freedv_get_bits_per_modem_frame       (struct freedv *f){return f->bits_per_modem_frame;}
int freedv_get_uncorrected_errors         (struct freedv *f) {return f->rx_status & RX_BIT_ERRORS;}

int freedv_get_n_max_speech_samples(struct freedv *f) {
    /* When "passing through" demod samples to the speech output
       (e.g. no sync and squelch off) f->nin bounces around with
       timing variations.  So is is possible we may return
       freedv_get_n_max_modem_samples() via the speech_output[]
       array */
    if (freedv_get_n_max_modem_samples(f) > f->n_speech_samples)
        return freedv_get_n_max_modem_samples(f);
    else
        return f->n_speech_samples;
}

// Now dummy obsolete call
int freedv_get_sync_interleaver(struct freedv *f) {
    return 1;
}

int freedv_get_sz_error_pattern(struct freedv *f) 
{
    if (FDV_MODE_ACTIVE( FREEDV_MODE_700C, f->mode)) {
        /* if diversity disabled callback sends error pattern for upper and lower carriers */
        return f->sz_error_pattern * (2 - f->test_frames_diversity);
    } else {
        return f->sz_error_pattern;
    }
}

// Get modem status, scatter/eye diagram for plotting, other goodies
void freedv_get_modem_extended_stats(struct freedv *f, struct MODEM_STATS *stats)
{
    if (FDV_MODE_ACTIVE( FREEDV_MODE_1600, f->mode))
        fdmdv_get_demod_stats(f->fdmdv, stats);

    if (FDV_MODE_ACTIVE( FREEDV_MODE_2400A, f->mode) || FDV_MODE_ACTIVE( FREEDV_MODE_800XA, f->mode)) {
        fsk_get_demod_stats(f->fsk, stats);
        float EbNodB = stats->snr_est;                       /* fsk demod actually estimates Eb/No     */
        stats->snr_est = EbNodB + 10.0f*log10f(800.0f/3000.0f); /* so convert to SNR Rb=800, noise B=3000 */
    }

    if (FDV_MODE_ACTIVE( FREEDV_MODE_2400B, f->mode)) {
        fmfsk_get_demod_stats(f->fmfsk, stats);
    }
    
    if (FDV_MODE_ACTIVE( FREEDV_MODE_700C, f->mode)) {
        cohpsk_get_demod_stats(f->cohpsk, stats);
    }
    
    if (FDV_MODE_ACTIVE( FREEDV_MODE_700D, f->mode) || FDV_MODE_ACTIVE( FREEDV_MODE_2020, f->mode)) {
        ofdm_get_demod_stats(f->ofdm, stats);
    }
}

