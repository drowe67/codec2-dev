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
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "codec2.h"
#include "codec2_fdmdv.h"
#include "fdmdv_internal.h"
#include "golay23.h"
#include "varicode.h"
#include "freedv_api.h"
#include "freedv_api_internal.h"
#include "comp_prim.h"

#define VERSION     11    /* The API version number.  The first version
                           is 10.  Increment if the API changes in a
                           way that would require changes by the API
                           user. */
/*
 * Version 10   Initial version August 2, 2015.
 * Version 11   September 2015
 *              Added: freedv_zero_total_bit_errors(), freedv_get_sync()
 *              Changed all input and output sample rates to 8000 sps.  Rates for FREEDV_MODE_700 and 700B were 7500.
 */

#define NORM_PWR  1.74   /* experimentally derived fudge factor so 1600 and
                            700 mode have the same tx power */

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

    if ((mode != FREEDV_MODE_1600) && (mode != FREEDV_MODE_700) && (mode != FREEDV_MODE_700B))
        return NULL;

    f = (struct freedv*)malloc(sizeof(struct freedv));
    if (f == NULL)
        return NULL;

    f->mode = mode;
    f->test_frames = f->smooth_symbols = 0;
    f->freedv_put_error_pattern = NULL;
    f->error_pattern_callback_state = NULL;

    if (mode == FREEDV_MODE_1600) {
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
        f->fdmdv_bits = (int*)malloc(nbit*sizeof(int));
        if (f->fdmdv_bits == NULL)
            return NULL;
        nbit = 2*fdmdv_bits_per_frame(f->fdmdv);
        f->tx_bits = (int*)malloc(nbit*sizeof(int));
        f->rx_bits = (int*)malloc(nbit*sizeof(int));
        if ((f->tx_bits == NULL) || (f->rx_bits == NULL))
            return NULL;
        f->evenframe = 0;
        f->sz_error_pattern = fdmdv_error_pattern_size(f->fdmdv);
    }

#ifndef CORTEX_M4
    if ((mode == FREEDV_MODE_700) || (mode == FREEDV_MODE_700B)) {
        f->snr_squelch_thresh = 0.0;
        f->squelch_en = 0;
        if (mode == FREEDV_MODE_700)
            codec2_mode = CODEC2_MODE_700;
        else
            codec2_mode = CODEC2_MODE_700B;
        f->cohpsk = cohpsk_create();
        f->nin = COHPSK_NOM_SAMPLES_PER_FRAME;
        f->n_nat_modem_samples = COHPSK_NOM_SAMPLES_PER_FRAME;          // native modem samples as used by the modem
        f->n_nom_modem_samples = f->n_nat_modem_samples * 8000 / 7500;  // number of samples after native samples are interpolated to 8000 sps
        f->n_max_modem_samples = COHPSK_MAX_SAMPLES_PER_FRAME * 8000 / 7500 + 1;
        f->modem_sample_rate = COHPSK_FS;                /* note wierd sample rate */
        f->clip = 1;
        nbit = COHPSK_BITS_PER_FRAME;
        f->tx_bits = (int*)malloc(nbit*sizeof(int));
        if (f->tx_bits == NULL)
            return NULL;
        f->sz_error_pattern = cohpsk_error_pattern_size();
    }
#endif

    f->test_frame_sync_state = 0;
    f->total_bits = 0;
    f->total_bit_errors = 0;

    f->codec2 = codec2_create(codec2_mode);
    if (f->codec2 == NULL)
        return NULL;
    if (mode == FREEDV_MODE_1600)
        f->n_speech_samples = codec2_samples_per_frame(f->codec2);
    if ((mode == FREEDV_MODE_700) || (mode == FREEDV_MODE_700B))
        f->n_speech_samples = 2*codec2_samples_per_frame(f->codec2);
    f->prev_rx_bits = (float*)malloc(sizeof(float)*2*codec2_bits_per_frame(f->codec2));
    if (f->prev_rx_bits == NULL)
        return NULL;

    nbit = codec2_bits_per_frame(f->codec2);
    nbyte = (nbit + 7) / 8;
    f->packed_codec_bits = (unsigned char*)malloc(nbyte*sizeof(char));
    if (mode == FREEDV_MODE_1600)
        f->codec_bits = (int*)malloc(nbit*sizeof(int));
    if ((mode == FREEDV_MODE_700) || (mode == FREEDV_MODE_700B))
        f->codec_bits = (int*)malloc(COHPSK_BITS_PER_FRAME*sizeof(int));

    if ((f->packed_codec_bits == NULL) || (f->codec_bits == NULL))
        return NULL;

    if ((mode == FREEDV_MODE_700) || (mode == FREEDV_MODE_700B)) {        // change modem rates to 8000 sps
        f->ptFilter7500to8000 = (struct quisk_cfFilter *)malloc(sizeof(struct quisk_cfFilter));
        f->ptFilter8000to7500 = (struct quisk_cfFilter *)malloc(sizeof(struct quisk_cfFilter));
        quisk_filt_cfInit(f->ptFilter8000to7500, quiskFilt120t480, sizeof(quiskFilt120t480)/sizeof(float));
        quisk_filt_cfInit(f->ptFilter7500to8000, quiskFilt120t480, sizeof(quiskFilt120t480)/sizeof(float));
    }
    else {
        f->ptFilter7500to8000 = NULL;
        f->ptFilter8000to7500 = NULL;
    }

    varicode_decode_init(&f->varicode_dec_states, 1);
    f->nvaricode_bits = 0;
    f->varicode_bit_index = 0;
    f->freedv_get_next_tx_char = NULL;
    f->freedv_put_next_rx_char = NULL;

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

    free(freedv->prev_rx_bits);
    free(freedv->packed_codec_bits);
    free(freedv->codec_bits);
    free(freedv->tx_bits);
    if (freedv->mode == FREEDV_MODE_1600)
        fdmdv_destroy(freedv->fdmdv);
#ifndef CORTEX_M4
    if (freedv->mode == FREEDV_MODE_700)
        cohpsk_destroy(freedv->cohpsk);
#endif
    codec2_destroy(freedv->codec2);
    if (freedv->ptFilter8000to7500) {
        quisk_filt_destroy(freedv->ptFilter8000to7500);
        free(freedv->ptFilter8000to7500);
        freedv->ptFilter8000to7500 = NULL;
    }
    if (freedv->ptFilter7500to8000) {
        quisk_filt_destroy(freedv->ptFilter7500to8000);
        free(freedv->ptFilter7500to8000);
        freedv->ptFilter7500to8000 = NULL;
    }
    free(freedv);
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

void freedv_tx(struct freedv *f, short mod_out[], short speech_in[]) {
    assert(f != NULL);
    COMP tx_fdm[f->n_nom_modem_samples];
    int  i;

    freedv_comptx(f, tx_fdm, speech_in);

    for(i=0; i<f->n_nom_modem_samples; i++)
        mod_out[i] = tx_fdm[i].real;
}

/* complex valued output, useful for suitable for single sided freq shifting */

void freedv_comptx(struct freedv *f, COMP mod_out[], short speech_in[]) {
    assert(f != NULL);
    int    bit, byte, i, j, k;
    int    bits_per_codec_frame, bits_per_modem_frame;
    int    data, codeword1, data_flag_index, nspare;
    COMP   tx_fdm[f->n_nat_modem_samples];

    assert((f->mode == FREEDV_MODE_1600) || (f->mode == FREEDV_MODE_700) || (f->mode == FREEDV_MODE_700B));

    if (f->mode == FREEDV_MODE_1600) {
        bits_per_codec_frame = codec2_bits_per_frame(f->codec2);
        bits_per_modem_frame = fdmdv_bits_per_frame(f->fdmdv);

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


#ifndef CORTEX_M4
    if ((f->mode == FREEDV_MODE_700) || (f->mode == FREEDV_MODE_700B)) {
        bits_per_codec_frame = codec2_bits_per_frame(f->codec2);
        bits_per_modem_frame = COHPSK_BITS_PER_FRAME;

        for (j=0; j<bits_per_modem_frame; j+=bits_per_codec_frame) {
            codec2_encode(f->codec2, f->packed_codec_bits, speech_in);
            speech_in += codec2_samples_per_frame(f->codec2);

            /* unpack bits, MSB first */

            bit = 7; byte = 0;
            for(i=0; i<bits_per_codec_frame; i++) {
                f->codec_bits[j+i] = (f->packed_codec_bits[byte] >> bit) & 0x1;
                bit--;
                if (bit < 0) {
                    bit = 7;
                    byte++;
                }
            }

            // spare bits in frame that codec defines.  Use these spare
            // bits/frame to send txt messages

            if (f->mode == FREEDV_MODE_700)
                nspare = 2;
            else
                nspare = 1; // Just one spare bit for FREEDV_MODE_700B

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

        cohpsk_mod(f->cohpsk, tx_fdm, f->codec_bits);
        if (f->clip)
            cohpsk_clip(tx_fdm);
        for(i=0; i<f->n_nat_modem_samples; i++)
            mod_out[i] = fcmult(FDMDV_SCALE*NORM_PWR, tx_fdm[i]);
        i = quisk_cfInterpDecim(mod_out, f->n_nat_modem_samples, f->ptFilter7500to8000, 16, 15);
        //assert(i == f->n_nom_modem_samples);
        // Caution: assert fails if f->n_nat_modem_samples * 16.0 / 15.0 is not an integer
    }
#endif
}

int freedv_nin(struct freedv *f) {
    if ((f->mode == FREEDV_MODE_700) || (f->mode == FREEDV_MODE_700B))
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

  1600 mode: When out of sync, it the number of output speech samples
  returned will be freedv_nin(). When in sync to a valid FreeDV 1600
  signal, the number of output speech samples will alternate between
  freedv_get_n_speech_samples() and 0.

  700 and 700B mode: The number of output speech samples returned will
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
    COMP rx_fdm[f->n_max_modem_samples];
    int i;
    int nin = freedv_nin(f);

    assert(nin <= f->n_max_modem_samples);

    for(i=0; i<nin; i++) {
        rx_fdm[i].real = (float)demod_in[i];
        rx_fdm[i].imag = 0.0;
    }

    return freedv_comprx(f, speech_out, rx_fdm);
}


// float input samples version

int freedv_floatrx(struct freedv *f, short speech_out[], float demod_in[]) {
    assert(f != NULL);
    COMP rx_fdm[f->n_max_modem_samples];
    int  i;
    int nin = freedv_nin(f);

    assert(nin <= f->n_max_modem_samples);

    for(i=0; i<nin; i++) {
        rx_fdm[i].real = demod_in[i];
        rx_fdm[i].imag = 0;
    }

    return freedv_comprx(f, speech_out, rx_fdm);
}

// complex input samples version

int freedv_comprx(struct freedv *f, short speech_out[], COMP demod_in[]) {
    assert(f != NULL);
    int                 bits_per_codec_frame, bytes_per_codec_frame, bits_per_fdmdv_frame;
    int                 i, j, bit, byte, nin_prev, nout, k;
    int                 recd_codeword, codeword1, data_flag_index, n_ascii, nspare;
    short               abit[1];
    char                ascii_out;
    int                 reliable_sync_bit;

    assert(f->nin <= f->n_max_modem_samples);

    bits_per_codec_frame  = codec2_bits_per_frame(f->codec2);
    bytes_per_codec_frame = (bits_per_codec_frame + 7) / 8;
    nout = f->n_speech_samples;

    if (f->mode == FREEDV_MODE_1600) {

        for(i=0; i<f->nin; i++)
            demod_in[i] = fcmult(1.0/FDMDV_SCALE, demod_in[i]);

        bits_per_fdmdv_frame  = fdmdv_bits_per_frame(f->fdmdv);

        nin_prev = f->nin;
        fdmdv_demod(f->fdmdv, f->fdmdv_bits, &reliable_sync_bit, demod_in, &f->nin);
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

                    codec2_decode(f->codec2, speech_out, f->packed_codec_bits);
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
                    for(i=0; i<f->n_speech_samples; i++)
                        speech_out[i] = 0;
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
                for(i=0; i<nin_prev; i++)
                    speech_out[i] = FDMDV_SCALE*demod_in[i].real;
            }
            else {
                for(i=0; i<nin_prev; i++)
                    speech_out[i] = 0;
             }
            //fprintf(stderr, "%d %d %d\n", nin_prev, speech_out[0], speech_out[nin_prev-1]);
            nout = nin_prev;
        }
    }


#ifndef CORTEX_M4
    if ((f->mode == FREEDV_MODE_700) || (f->mode == FREEDV_MODE_700B)) {
        float rx_bits[COHPSK_BITS_PER_FRAME];
        int   sync;

        // echo samples back out as default (say if sync not found)

        for(i=0; i<freedv_nin(f); i++)
            speech_out[i] = demod_in[i].real;

        i = quisk_cfInterpDecim(demod_in, freedv_nin(f), f->ptFilter8000to7500, 15, 16);
        //if (i != f->nin)
        //    printf("freedv_comprx decimation: input %d output %d\n", freedv_nin(f), i);

        for(i=0; i<f->nin; i++)
            demod_in[i] = fcmult(1.0/FDMDV_SCALE, demod_in[i]);

        nin_prev = f->nin;
        cohpsk_demod(f->cohpsk, rx_bits, &sync, demod_in, &f->nin);
        f->sync = sync;
        cohpsk_get_demod_stats(f->cohpsk, &f->stats);
        f->snr_est = f->stats.snr_est;

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

                for (j=0; j<COHPSK_BITS_PER_FRAME; j+=bits_per_codec_frame) {

                    /* extract txt msg data bit(s) */

                    if (f->mode == FREEDV_MODE_700)
                        nspare = 2;
                    else
                        nspare = 1;

                    for(k=0; k<nspare; k++)  {
                        abit[0] = rx_bits[data_flag_index+j+k] < 0.0;

                        n_ascii = varicode_decode(&f->varicode_dec_states, &ascii_out, abit, 1, 1);
                        if (n_ascii && (f->freedv_put_next_rx_char != NULL)) {
                            (*f->freedv_put_next_rx_char)(f->callback_state, ascii_out);
                        }
                    }

                    /* pack bits, MSB received first */

                    bit = 7; byte = 0;
                    memset(f->packed_codec_bits, 0, bytes_per_codec_frame);
                    for(i=0; i<bits_per_codec_frame; i++) {
                        f->packed_codec_bits[byte] |= ((rx_bits[j+i] < 0.0) << bit);
                        bit--;
                        if (bit < 0) {
                            bit = 7;
                            byte++;
                        }
                    }

                    codec2_decode(f->codec2, speech_out, f->packed_codec_bits);

                    if (f->squelch_en && (f->stats.snr_est < f->snr_squelch_thresh)) {
                        for(i=0; i<f->n_speech_samples; i++)
                            speech_out[i] = 0;
                    }
                    speech_out += codec2_samples_per_frame(f->codec2);
                }
                nout = f->n_speech_samples;
            }
            else {
                short error_pattern[COHPSK_BITS_PER_FRAME];
                int   bit_errors;

                /* test data, lets see if we can sync to the test data sequence */

                cohpsk_put_test_bits(f->cohpsk, &f->test_frame_sync_state, error_pattern, &bit_errors, rx_bits);
                if (f->test_frame_sync_state) {
                    f->total_bit_errors += bit_errors;
                    f->total_bits       += COHPSK_BITS_PER_FRAME;
                    if (f->freedv_put_error_pattern != NULL) {
                        (*f->freedv_put_error_pattern)(f->error_pattern_callback_state, error_pattern, COHPSK_BITS_PER_FRAME);
                    }
                }

                for(i=0; i<f->n_speech_samples; i++)
                    speech_out[i] = 0;
                nout = f->n_speech_samples;
            }

        }


        /* no valid FreeDV signal - squelch output */

        if (sync == 0) {
            nout = freedv_nin(f);
            if (f->squelch_en) {
                for(i=0; i<nout; i++)
                    speech_out[i] = 0;
            }
        }

    }
#endif

    //fprintf(stderr,"freedv_nin(f): %d nout: %d\n", freedv_nin(f), nout);
    return nout;
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
    f->freedv_put_next_rx_char = rx;
    f->freedv_get_next_tx_char = tx;
    f->callback_state = state;
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
    if (f->mode == FREEDV_MODE_1600)
        fdmdv_get_demod_stats(f->fdmdv, &f->stats);
#ifndef CORTEX_M4
    if ((f->mode == FREEDV_MODE_700) || (f->mode == FREEDV_MODE_700B))
        cohpsk_get_demod_stats(f->cohpsk, &f->stats);
#endif
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
void freedv_set_squelch_en                (struct freedv *f, int val) {f->squelch_en = val;}
void freedv_set_total_bit_errors          (struct freedv *f, int val) {f->total_bit_errors = val;}
void freedv_set_total_bits                (struct freedv *f, int val) {f->total_bits = val;}
void freedv_set_clip                       (struct freedv *f, int val) {f->clip = val;}
void freedv_set_varicode_code_num         (struct freedv *f, int val) {varicode_set_code_num(&f->varicode_dec_states, val);}

// Set floats
void freedv_set_snr_squelch_thresh        (struct freedv *f, float val) {f->snr_squelch_thresh = val;}

void freedv_set_callback_error_pattern    (struct freedv *f, freedv_calback_error_pattern cb, void *state)
{
    f->freedv_put_error_pattern = cb;
    f->error_pattern_callback_state = state;
}

/*---------------------------------------------------------------------------*\

  FUNCTIONS...: freedv_get_*
  AUTHOR......: Jim Ahlstrom
  DATE CREATED: 28 July 2015

  Get some parameters from FreeDV.  It is possible to write a macro using ## for
  this, but I wasn't sure it would be 100% portable.

\*---------------------------------------------------------------------------*/

// Get integers
int freedv_get_mode                       (struct freedv *f) {return f->mode;}
int freedv_get_test_frames                (struct freedv *f) {return f->test_frames;}
int freedv_get_n_speech_samples           (struct freedv *f) {return f->n_speech_samples;}
int freedv_get_modem_sample_rate          (struct freedv *f) {return FS;}
int freedv_get_n_max_modem_samples        (struct freedv *f) {return f->n_max_modem_samples;}
int freedv_get_n_nom_modem_samples        (struct freedv *f) {return f->n_nom_modem_samples;}
int freedv_get_total_bits                 (struct freedv *f) {return f->total_bits;}
int freedv_get_total_bit_errors           (struct freedv *f) {return f->total_bit_errors;}
int freedv_get_sync                       (struct freedv *f) {return  f->stats.sync;}
int freedv_get_sz_error_pattern           (struct freedv *f) {return  f->sz_error_pattern;}
// Get floats

struct CODEC2 *freedv_get_codec2	(struct freedv *f){return  f->codec2;}

void freedv_get_modem_extended_stats(struct freedv *f, struct MODEM_STATS *stats)
{
    if (f->mode == FREEDV_MODE_1600)
        fdmdv_get_demod_stats(f->fdmdv, stats);
#ifndef CORTEX_M4
    if ((f->mode == FREEDV_MODE_700) || (f->mode == FREEDV_MODE_700B))
        cohpsk_get_demod_stats(f->cohpsk, stats);
#endif
}

/*--  Functions below this line are private, and not meant for public use  --*/
/*---------------------------------------------------------------------------*\

  FUNCTIONS...: quisk_filt_cfInit
  AUTHOR......: Jim Ahlstrom
  DATE CREATED: 27 August 2015

  Initialize a FIR filter that will be used to change sample rates.  These rate
  changing filters were copied from Quisk and modified for float samples.

\*---------------------------------------------------------------------------*/

static void quisk_filt_cfInit(struct quisk_cfFilter * filter, float * coefs, int taps)
{    // Prepare a new filter using coefs and taps.  Samples are complex.
    filter->dCoefs = coefs;
    filter->cSamples = (COMP *)malloc(taps * sizeof(COMP));
    memset(filter->cSamples, 0, taps * sizeof(COMP));
    filter->ptcSamp = filter->cSamples;
    filter->nTaps = taps;
    filter->cBuf = NULL;
    filter->nBuf = 0;
    filter->decim_index = 0;
}

/*---------------------------------------------------------------------------*\

  FUNCTIONS...: quisk_filt_destroy
  AUTHOR......: Jim Ahlstrom
  DATE CREATED: 27 August 2015

  Destroy the FIR filter and free all resources.

\*---------------------------------------------------------------------------*/

static void quisk_filt_destroy(struct quisk_cfFilter * filter)
{
    if (filter->cSamples) {
        free(filter->cSamples);
        filter->cSamples = NULL;
    }
    if (filter->cBuf) {
        free(filter->cBuf);
        filter->cBuf = NULL;
    }
}

/*---------------------------------------------------------------------------*\

  FUNCTIONS...: quisk_cfInterpDecim
  AUTHOR......: Jim Ahlstrom
  DATE CREATED: 27 August 2015

  Take an array of samples cSamples of length count, multiply the sample rate
  by interp, and then divide the sample rate by decim.  Return the new number
  of samples.  Each specific interp and decim will require its own custom
  FIR filter.

\*---------------------------------------------------------------------------*/

static int quisk_cfInterpDecim(COMP * cSamples, int count, struct quisk_cfFilter * filter, int interp, int decim)
{   // Interpolate by interp, and then decimate by decim.
    // This uses the float coefficients of filter (not the complex).  Samples are complex.
    int i, k, nOut;
    float * ptCoef;
    COMP * ptSample;
    COMP csample;

    if (count > filter->nBuf) {    // increase size of sample buffer
        filter->nBuf = count * 2;
        if (filter->cBuf)
            free(filter->cBuf);
        filter->cBuf = (COMP *)malloc(filter->nBuf * sizeof(COMP));
    }
    memcpy(filter->cBuf, cSamples, count * sizeof(COMP));
    nOut = 0;
    for (i = 0; i < count; i++) {
        // Put samples into buffer left to right.  Use samples right to left.
        *filter->ptcSamp = filter->cBuf[i];
        while (filter->decim_index < interp) {
            ptSample = filter->ptcSamp;
            ptCoef = filter->dCoefs + filter->decim_index;
            csample.real = 0;
            csample.imag = 0;
            for (k = 0; k < filter->nTaps / interp; k++, ptCoef += interp) {
                csample.real += (*ptSample).real * *ptCoef;
                csample.imag += (*ptSample).imag * *ptCoef;
                if (--ptSample < filter->cSamples)
                    ptSample = filter->cSamples + filter->nTaps - 1;
            }
            cSamples[nOut].real = csample.real * interp;
            cSamples[nOut].imag = csample.imag * interp;
            nOut++;
            filter->decim_index += decim;
        }
        if (++filter->ptcSamp >= filter->cSamples + filter->nTaps)
            filter->ptcSamp = filter->cSamples;
        filter->decim_index = filter->decim_index - interp;
    }
    return nOut;
}

