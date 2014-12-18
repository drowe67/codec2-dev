/*---------------------------------------------------------------------------*\
                                                                             
  FILE........: freedv_api.c
  AUTHOR......: David Rowe
  DATE CREATED: August 2014
                                                                             
  Library of API functions that implement FreeDV "modes", useful for
  embedding FreeDV in other programs.
      
  TODO:
    [X] speex tx/rx works
    [X] txt messages
    [ ] optional test tx framemode
                                                                       
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
#include <stdio.h>
#include <string.h>

#include "codec2.h"
#include "codec2_fdmdv.h"
#include "golay23.h"
#include "varicode.h"
#include "freedv_api.h"

/*---------------------------------------------------------------------------*\
                                                       
  FUNCTION....: freedv_open
  AUTHOR......: David Rowe			      
  DATE CREATED: 3 August 2014

  Call this first to initialise.  Returns NULL if initialisation fails
  (e.g. out of memory or mode not supported).

\*---------------------------------------------------------------------------*/

struct freedv *freedv_open(int mode) {
    struct freedv *f;
    int            Nc, codec2_mode, nbit, nbyte;
    
    if (mode != FREEDV_MODE_1600)
        return NULL;

    f = (struct freedv*)malloc(sizeof(struct freedv));
    if (f == NULL)
        return NULL;

    f->mode = mode;
    f->tx_sync_bit = 0;
    f->snr_thresh = 2.0;

    if (mode == FREEDV_MODE_1600) {
        Nc = 16;
        codec2_mode = CODEC2_MODE_1300;
    }

    f->codec2 = codec2_create(codec2_mode);
    if (f->codec2 == NULL)
        return NULL;

    f->fdmdv = fdmdv_create(Nc);
    if (f->fdmdv == NULL)
        return NULL;

    nbit = codec2_bits_per_frame(f->codec2);
    nbyte = (nbit + 7) / 8;
    f->packed_codec_bits = (unsigned char*)malloc(nbyte*sizeof(char));
    f->codec_bits = (int*)malloc(nbit*sizeof(int));

    nbit = 2*fdmdv_bits_per_frame(f->fdmdv);
    f->tx_bits = (int*)malloc(nbit*sizeof(int));
    f->rx_bits = (int*)malloc(nbit*sizeof(int));
    
    nbit = fdmdv_bits_per_frame(f->fdmdv);
    f->fdmdv_bits = (int*)malloc(nbit*sizeof(int));

    if ((f->packed_codec_bits == NULL) || (f->codec_bits == NULL) 
        || (f->tx_bits == NULL) || (f->rx_bits == NULL) || (f->fdmdv_bits == NULL))
        return NULL;

    varicode_decode_init(&f->varicode_dec_states, 1);
    f->nvaricode_bits = 0;
    f->varicode_bit_index = 0;
    f->freedv_get_next_tx_char = NULL;
    f->freedv_put_next_rx_char = NULL;

    golay23_init();
    f->total_bit_errors = 0;

    f->nin = FDMDV_NOM_SAMPLES_PER_FRAME;

    return f;
}

/*---------------------------------------------------------------------------*\
                                                       
  FUNCTION....: freedv_close
  AUTHOR......: David Rowe			      
  DATE CREATED: 3 August 2014

  Frees up memory.

\*---------------------------------------------------------------------------*/

void freedv_close(struct freedv *freedv) {
    free(freedv->packed_codec_bits);
    free(freedv->codec_bits);
    free(freedv->tx_bits);
    fdmdv_destroy(freedv->fdmdv);
    codec2_destroy(freedv->codec2);
    free(freedv);
}

/*---------------------------------------------------------------------------*\
                                                       
  FUNCTION....: freedv_tx
  AUTHOR......: David Rowe			      
  DATE CREATED: 3 August 2014

  Takes a frame of input speech samples, encodes and modulates them to produce
  a frame of modem samples that can be sent to the transmitter.

  speech_in[] and mod_out[] are sampled at 8 kHz, 16 bit shorts, and
  are always FREEDV_NSAMPLES long.  The speech_in[] level should be
  such that the peak speech level is between +/16384 and +/- 32767.
  mod_out[] will be scaled such that the peak level is just less than
  +/-32767.

  The FDM modem signal mod_out[] has a high crest factor (around
  12dB), however the energy and duration of the peaks is small.
  FreeDV is generally operated at a "backoff" of 6-8dB.  Adjust the
  power amplifier drive so that the average power is 6-8dB less than
  the peak power of the PA.  For example, on a radio rated at 100W PEP
  for SSB, the average FreeDV power is typically 20-25W.

\*---------------------------------------------------------------------------*/

void freedv_tx(struct freedv *f, short mod_out[], short speech_in[]) {
    int    bit, byte, i, j;
    int    bits_per_codec_frame, bits_per_fdmdv_frame;
    int    data, codeword1, data_flag_index;
    COMP   tx_fdm[2*FDMDV_NOM_SAMPLES_PER_FRAME];
     
    bits_per_codec_frame = codec2_bits_per_frame(f->codec2);
    bits_per_fdmdv_frame = fdmdv_bits_per_frame(f->fdmdv);

    codec2_encode(f->codec2, f->packed_codec_bits, speech_in);

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

    if (f->mode == FREEDV_MODE_1600) {
            
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

        //for(i=0; i<bits_per_codec_frame+12; i++)
        //    printf("%d\n", f->tx_bits[i]);
    }

    /* modulate even and odd frames */

    fdmdv_mod(f->fdmdv, tx_fdm, f->tx_bits, &f->tx_sync_bit);
    assert(f->tx_sync_bit == 1);

    fdmdv_mod(f->fdmdv, &tx_fdm[FDMDV_NOM_SAMPLES_PER_FRAME], &f->tx_bits[bits_per_fdmdv_frame], &f->tx_sync_bit);
    assert(f->tx_sync_bit == 0);

    for(i=0; i<2*FDMDV_NOM_SAMPLES_PER_FRAME; i++)
        mod_out[i] = FDMDV_SCALE * tx_fdm[i].real;

    assert(2*FDMDV_NOM_SAMPLES_PER_FRAME == FREEDV_NSAMPLES);
}

int freedv_nin(struct freedv *f) {
    return f->nin;
}

/*---------------------------------------------------------------------------*\
                                                       
  FUNCTION....: freedv_rx
  AUTHOR......: David Rowe			      
  DATE CREATED: 3 August 2014

  Takes a frame of samples from the radio receiver, demodulates them,
  then decodes them, producing a frame of decoded speech samples.  

  Both demod_in[] and speech_out[] are 16 bit shorts sampled at 8 kHz.
  The peak level of demod_in[] is not critical, as the demod works
  well over a wide range of amplitude scaling.  However it is best to
  avoid clipping (overload, or samples pinned to +/- 32767).  Suggest
  scaling so the peak (modem signal plus noise) is between +/-16384
  and +/-32767.  speech_out[] will peak at just less than +/-32767.

  To account for difference in the transmit and receive sample clock
  frequencies, the number of demod_in[] samples is time varying.  It
  is the responsibility of the caller to supply the correct number of
  samples.  Call freedv_nin() before each call to freedv_rx() to
  determine how many samples to pass to this function (see example).

  Returns the number of output speech samples available in
  speech_out[]. When in sync this will typically alternate between 0
  and FREEDV_NSAMPLES.  When out of sync, this will be f->nin.  

  When out of sync, this function echoes the demod_in[] samples to
  speech_out[].  This allows the user to listen to the channel, which
  is useful for tuning FreeDV signals or reception of non-FreeDV
  signals.

\*---------------------------------------------------------------------------*/

// short version

int freedv_rx(struct freedv *f, short speech_out[], short demod_in[]) {

float rxdata[FDMDV_MAX_SAMPLES_PER_FRAME];
int i;

    for(i=0; i<f->nin; i++)  rxdata[i] = (float)demod_in[i]/FDMDV_SCALE;

    return freedv_floatrx(f, speech_out, rxdata);
    
}


// float version

int freedv_floatrx(struct freedv *f, short speech_out[], float demod_in[]) {
    COMP                rx_fdm[FDMDV_MAX_SAMPLES_PER_FRAME];
    int                 bits_per_codec_frame, bytes_per_codec_frame, bits_per_fdmdv_frame;
    int                 reliable_sync_bit, i, j, bit, byte, nin_prev, nout;
    int                 recd_codeword, codeword1, data_flag_index, n_ascii;
    short               abit[1];
    char                ascii_out;

    bits_per_codec_frame  = codec2_bits_per_frame(f->codec2);
    bytes_per_codec_frame = (bits_per_codec_frame + 7) / 8;
    bits_per_fdmdv_frame  = fdmdv_bits_per_frame(f->fdmdv);

    for(i=0; i<f->nin; i++) {
        rx_fdm[i].real = demod_in[i];
        rx_fdm[i].imag = 0;
    }

    nin_prev = f->nin;
    fdmdv_demod(f->fdmdv, f->fdmdv_bits, &reliable_sync_bit, rx_fdm, &f->nin);
    fdmdv_get_demod_stats(f->fdmdv, &f->fdmdv_stats);

    if (f->fdmdv_stats.sync) {
        if (reliable_sync_bit == 0) {
            memcpy(f->rx_bits, f->fdmdv_bits, bits_per_fdmdv_frame*sizeof(int));
            nout = 0;
        }
        else {
            memcpy(&f->rx_bits[bits_per_fdmdv_frame], f->fdmdv_bits, bits_per_fdmdv_frame*sizeof(int));
   
            if (f->mode == FREEDV_MODE_1600) {
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

            codec2_decode(f->codec2, speech_out, f->packed_codec_bits);

            /* squelch if beneath SNR threshold */

            if (f->fdmdv_stats.snr_est < f->snr_thresh) {
                for(i=0; i<FREEDV_NSAMPLES; i++)
                    speech_out[i] = 0;
            }

            nout = FREEDV_NSAMPLES;
        }
    } /* if (sync) .... */
    else {
        /* if not in sync pass through analog samples */
        /* this lets us "hear" whats going on, e.g. during tuning */
        for(i=0; i<nin_prev; i++)
            speech_out[i] = demod_in[i];
        nout = nin_prev;
    }

    return nout;
}

