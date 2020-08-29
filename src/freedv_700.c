/*---------------------------------------------------------------------------*\

  FILE........: freedv_700.c
  AUTHOR......: David Rowe
  DATE CREATED: May 2020

  Functions that implement the various FreeDV 700 modes.

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
#include "codec2_fdmdv.h"
#include "varicode.h"
#include "freedv_api.h"
#include "freedv_api_internal.h"
#include "comp_prim.h"

#include "codec2_ofdm.h"
#include "ofdm_internal.h"
#include "mpdecode_core.h"
#include "gp_interleaver.h"
#include "ldpc_codes.h"
#include "interldpc.h"
#include "debug_alloc.h"
#include "filter.h"

extern char *ofdm_statemode[];

void freedv_700c_open(struct freedv *f) {
    f->snr_squelch_thresh = 0.0;
    f->squelch_en = 0;

    f->cohpsk = cohpsk_create();
    f->nin = f->nin_prev = COHPSK_NOM_SAMPLES_PER_FRAME;
    f->n_nat_modem_samples = COHPSK_NOM_SAMPLES_PER_FRAME;                       // native modem samples as used by the modem
    f->n_nom_modem_samples = f->n_nat_modem_samples * FREEDV_FS_8000 / COHPSK_FS;// number of samples after native samples are interpolated to 8000 sps
    f->n_max_modem_samples = COHPSK_MAX_SAMPLES_PER_FRAME * FREEDV_FS_8000 / COHPSK_FS + 1;
    f->modem_sample_rate = FREEDV_FS_8000;                                       // note weird sample rate tamed by resampling
    f->clip = 1;
    f->sz_error_pattern = cohpsk_error_pattern_size();
    f->test_frames_diversity = 1;

    f->ptFilter7500to8000 = (struct quisk_cfFilter *)MALLOC(sizeof(struct quisk_cfFilter));
    f->ptFilter8000to7500 = (struct quisk_cfFilter *)MALLOC(sizeof(struct quisk_cfFilter));
    quisk_filt_cfInit(f->ptFilter8000to7500, quiskFilt120t480, sizeof(quiskFilt120t480)/sizeof(float));
    quisk_filt_cfInit(f->ptFilter7500to8000, quiskFilt120t480, sizeof(quiskFilt120t480)/sizeof(float));

    f->speech_sample_rate = FREEDV_FS_8000;
    f->codec2 = codec2_create(CODEC2_MODE_700C); assert(f->codec2 != NULL);

    f->n_codec_frames = 2;
    f->n_speech_samples = f->n_codec_frames*codec2_samples_per_frame(f->codec2);
    f->bits_per_codec_frame = codec2_bits_per_frame(f->codec2);
    f->bits_per_modem_frame = f->n_codec_frames*codec2_bits_per_frame(f->codec2);
    assert(f->bits_per_modem_frame == COHPSK_BITS_PER_FRAME);

    f->tx_payload_bits = (uint8_t*)MALLOC(f->bits_per_modem_frame*sizeof(char)); assert(f->tx_payload_bits != NULL);
    f->rx_payload_bits = (uint8_t*)MALLOC(f->bits_per_modem_frame*sizeof(char)); assert(f->rx_payload_bits != NULL);
}

void freedv_comptx_700c(struct freedv *f, COMP mod_out[]) {
    int    i;
#ifdef NO_C99
    COMP   *tx_fdm = alloca(f->n_nat_modem_samples*sizeof(COMP));
#else
    COMP   tx_fdm[f->n_nat_modem_samples];
#endif
    int    tx_bits[COHPSK_BITS_PER_FRAME];

    /* earlier modems used one bit per int for unpacked bits */
    for(i=0; i<COHPSK_BITS_PER_FRAME; i++) tx_bits[i] = f->tx_payload_bits[i];
    
    /* optionally overwrite the codec bits with test frames */
    if (f->test_frames) {
        cohpsk_get_test_bits(f->cohpsk, tx_bits);
    }

    /* cohpsk modulator */
    cohpsk_mod(f->cohpsk, tx_fdm, tx_bits, COHPSK_BITS_PER_FRAME);
    if (f->clip) 
        cohpsk_clip(tx_fdm, COHPSK_CLIP, COHPSK_NOM_SAMPLES_PER_FRAME);
    for(i=0; i<f->n_nat_modem_samples; i++)
        mod_out[i] = fcmult(FDMDV_SCALE*NORM_PWR_COHPSK, tx_fdm[i]);
    i = quisk_cfInterpDecim((complexf_t *)mod_out, f->n_nat_modem_samples, f->ptFilter7500to8000, 16, 15);
}

void freedv_700d_open(struct freedv *f) {
    f->snr_squelch_thresh = 0.0;
    f->squelch_en = 0;

    f->ofdm = ofdm_create(NULL);
    assert(f->ofdm != NULL);
        
    struct OFDM_CONFIG *ofdm_config = ofdm_get_config_param(f->ofdm);
    f->ofdm_bitsperframe = ofdm_get_bits_per_frame(f->ofdm);
    f->ofdm_nuwbits = (ofdm_config->ns - 1) * ofdm_config->bps - ofdm_config->txtbits;
    f->ofdm_ntxtbits = ofdm_config->txtbits;

    f->ldpc = (struct LDPC*)MALLOC(sizeof(struct LDPC));
    assert(f->ldpc != NULL);

    ldpc_codes_setup(f->ldpc, "HRA_112_112");
#ifdef __EMBEDDED__
    f->ldpc->max_iter = 10; /* limit LDPC decoder iterations to limit CPU load */
#endif
    /* Code length 224 divided by 2 bits per symbol = 112 symbols per frame */
    int coded_syms_per_frame = f->ldpc->coded_bits_per_frame/f->ofdm->bps;
        
    f->modem_frame_count_tx = f->modem_frame_count_rx = 0;
        
    f->codeword_symbols = (COMP*)MALLOC(sizeof(COMP) * coded_syms_per_frame);
    assert(f->codeword_symbols != NULL);

    f->codeword_amps = (float*)MALLOC(sizeof(float) * coded_syms_per_frame);
    assert(f->codeword_amps != NULL);

    f->nin = f->nin_prev = ofdm_get_samples_per_frame(f->ofdm);
    f->n_nat_modem_samples = ofdm_get_samples_per_frame(f->ofdm);
    f->n_nom_modem_samples = ofdm_get_samples_per_frame(f->ofdm);
    f->n_max_modem_samples = ofdm_get_max_samples_per_frame(f->ofdm);
    f->modem_sample_rate = f->ofdm->config.fs;
    f->clip = 0;
    f->sz_error_pattern = f->ofdm_bitsperframe;

    f->tx_bits = NULL; /* not used for 700D */

#ifndef __EMBEDDED__
    /* tx BPF off on embedded platforms, as it consumes significant CPU */
    ofdm_set_tx_bpf(f->ofdm, 1);
#endif
        
    f->speech_sample_rate = FREEDV_FS_8000;
    f->codec2 = codec2_create(CODEC2_MODE_700C); assert(f->codec2 != NULL);
    /* should be exactly an integer number of Codec 2 frames in a OFDM modem frame */
    assert((f->ldpc->data_bits_per_frame % codec2_bits_per_frame(f->codec2)) == 0);

    f->n_codec_frames = f->ldpc->data_bits_per_frame/codec2_bits_per_frame(f->codec2);
    f->n_speech_samples = f->n_codec_frames*codec2_samples_per_frame(f->codec2);
    f->bits_per_codec_frame = codec2_bits_per_frame(f->codec2);
    f->bits_per_modem_frame = f->n_codec_frames*f->bits_per_codec_frame;

    f->tx_payload_bits = (unsigned char*)MALLOC(f->bits_per_modem_frame);
    assert(f->tx_payload_bits != NULL);
    f->rx_payload_bits = (unsigned char*)MALLOC(f->bits_per_modem_frame);
    assert(f->rx_payload_bits != NULL);
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

void freedv_comptx_700d(struct freedv *f, COMP mod_out[]) {
    int    i, k;
    int    nspare;
 
    // Generate Varicode txt bits. Txt bits in OFDM frame come just
    // after Unique Word (UW).  Txt bits aren't protected by FEC.

    nspare = f->ofdm_ntxtbits;
#ifdef NO_C99
    uint8_t *txt_bits = alloca(nspare*sizeof(uint8_t));
#else
    uint8_t txt_bits[nspare];
#endif

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
        else txt_bits[k] = 0;
    }

    /* optionally replace codec payload bits with test frames known to rx */

    if (f->test_frames) {
#ifdef NO_C99
        uint8_t *payload_data_bits = alloca(f->bits_per_modem_frame);
#else
        uint8_t payload_data_bits[f->bits_per_modem_frame];
#endif
        ofdm_generate_payload_data_bits(payload_data_bits, f->bits_per_modem_frame);

        for (i = 0; i < f->bits_per_modem_frame; i++) {
            f->tx_payload_bits[i] = payload_data_bits[i];
        }
    }

    /* OK now ready to LDPC encode, interleave, and OFDM modulate */
    
#ifdef NO_C99
    complexf_t *tx_sams = alloca(f->n_nat_modem_samples*sizeof(complexf_t));
#else
    complexf_t tx_sams[f->n_nat_modem_samples];
#endif
    COMP asam;
    
    ofdm_ldpc_interleave_tx(f->ofdm, f->ldpc, tx_sams, f->tx_payload_bits, txt_bits);

    for(i=0; i< f->n_nat_modem_samples; i++) {
        asam.real = crealf(tx_sams[i]);
        asam.imag = cimagf(tx_sams[i]);
        mod_out[i] = fcmult(OFDM_AMP_SCALE * NORM_PWR_OFDM, asam);
    }

    if (f->clip) {
        cohpsk_clip(mod_out, OFDM_CLIP, f->n_nat_modem_samples);
    }
}

int freedv_comprx_700c(struct freedv *f, COMP demod_in_8kHz[]) {
    int   i;
    int   sync;

    int rx_status = 0;
    
    // quisk_cfInterpDecim() modifies input data so lets make a copy just in case there
    // is no sync and we need to echo inpout to output

    // freedv_nin(f): input samples at Fs=8000 Hz
    // f->nin: input samples at Fs=7500 Hz
    
#ifdef NO_C99
    COMP *demod_in = alloca(freedv_nin(f)*sizeof(COMP));
#else
    COMP demod_in[freedv_nin(f)];
#endif
    
    for(i=0; i<freedv_nin(f); i++)
        demod_in[i] = demod_in_8kHz[i];

    i = quisk_cfInterpDecim((complexf_t *)demod_in, freedv_nin(f), f->ptFilter8000to7500, 15, 16);

    for(i=0; i<f->nin; i++)
        demod_in[i] = fcmult(1.0/FDMDV_SCALE, demod_in[i]);
    
    float rx_soft_bits[COHPSK_BITS_PER_FRAME];
    
    cohpsk_demod(f->cohpsk, rx_soft_bits, &sync, demod_in, &f->nin);
    
    for(i=0; i<f->bits_per_modem_frame; i++)
        f->rx_payload_bits[i] = rx_soft_bits[i] < 0.0f;
    
    f->sync = sync;
    cohpsk_get_demod_stats(f->cohpsk, &f->stats);
    f->snr_est = f->stats.snr_est;

    if (sync) {
        rx_status = RX_SYNC;
        
        if (f->test_frames == 0) {
            rx_status |= RX_BITS;
        }
        else {

            if (f->test_frames_diversity) {
                /* normal operation - error pattern on frame after diveristy combination */
                short error_pattern[COHPSK_BITS_PER_FRAME];
                int   bit_errors;

                /* test data, lets see if we can sync to the test data sequence */

                char rx_bits_char[COHPSK_BITS_PER_FRAME];
                for(i=0; i<COHPSK_BITS_PER_FRAME; i++)
                    rx_bits_char[i] = rx_soft_bits[i] < 0.0;
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

                /* combine total errors and call callback */

                if (f->test_frame_sync_state && f->test_frame_sync_state_upper) {
                    f->total_bit_errors += bit_errors_lower + bit_errors_upper;
                    f->total_bits       += 2*COHPSK_BITS_PER_FRAME;
                    if (f->freedv_put_error_pattern != NULL) {
                        (*f->freedv_put_error_pattern)(f->error_pattern_callback_state, error_pattern, 2*COHPSK_BITS_PER_FRAME);
                    }
                }

            }
        }

    }

    return rx_status;
}

/*
  700D demod function that can support complex (float) or real (short)
  samples.  The real short samples are useful for low memory overhead,
  such at the SM1000.
*/

int freedv_comp_short_rx_700d(struct freedv *f, void *demod_in_8kHz, int demod_in_is_short, float gain) {
    int   i, k;
    int   n_ascii;
    char  ascii_out;
    struct OFDM *ofdm = f->ofdm;
    struct LDPC *ldpc = f->ldpc;
    
    int    data_bits_per_frame = ldpc->data_bits_per_frame;
    int    coded_bits_per_frame = ldpc->coded_bits_per_frame;
    int    coded_syms_per_frame = ldpc->coded_bits_per_frame/ofdm->bps;
    COMP  *codeword_symbols = f->codeword_symbols;
    float *codeword_amps = f->codeword_amps;
#ifdef NO_C99
    int    *rx_bits = alloca(f->ofdm_bitsperframe*sizeof(int));
    short  *txt_bits = alloca(f->ofdm_ntxtbits*sizeof(short));
    COMP   *payload_syms = alloca(coded_syms_per_frame*sizeof(COMP));
    float  *payload_amps = alloca(coded_syms_per_frame*sizeof(float));
    uint8_t *rx_uw = alloca(f->ofdm_nuwbits);
#else
    int    rx_bits[f->ofdm_bitsperframe];
    short  txt_bits[f->ofdm_ntxtbits];
    COMP   payload_syms[coded_syms_per_frame];
    float  payload_amps[coded_syms_per_frame];
    uint8_t rx_uw[f->ofdm_nuwbits];
#endif
   
    int    Nerrs_raw = 0;
    int    Nerrs_coded = 0;
    int    iter = 0;
    int    parityCheckCount = 0;

    float new_gain = gain / OFDM_AMP_SCALE;
    
    assert((demod_in_is_short == 0) || (demod_in_is_short == 1));

    f->sync = f->stats.sync = 0;
    int rx_status = 0;
    
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
        rx_status |= RX_SYNC;
        if (ofdm->sync_state == trial) rx_status |= RX_TRIAL_SYNC;
    
        if (demod_in_is_short)
            ofdm_demod_shorts(ofdm, rx_bits, (short*)demod_in_8kHz, new_gain);
        else
            ofdm_demod(ofdm, rx_bits, (COMP*)demod_in_8kHz);

        ofdm_extract_uw(ofdm, ofdm->rx_np, ofdm->rx_amp, rx_uw);
        ofdm_disassemble_qpsk_modem_packet(ofdm, ofdm->rx_np, ofdm->rx_amp, payload_syms, payload_amps, txt_bits);

        f->sync = 1;
        ofdm_get_demod_stats(f->ofdm, &f->stats);
        f->snr_est = f->stats.snr_est;

        assert((f->ofdm_nuwbits+f->ofdm_ntxtbits+coded_bits_per_frame) == f->ofdm_bitsperframe);

        /* newest symbols at end of buffer (uses final i from last loop), note we 
           change COMP formats from what modem uses internally */
                
        for(i=0; i< coded_syms_per_frame; i++) {
            codeword_symbols[i] = payload_syms[i];
            codeword_amps[i]    = payload_amps[i];
         }
               
        /* run de-interleaver */
                
#ifdef NO_C99
        COMP  *codeword_symbols_de = alloca(coded_syms_per_frame*sizeof(COMP));
        float *codeword_amps_de = alloca(coded_syms_per_frame*sizeof(float));
        float *llr = alloca(coded_bits_per_frame*sizeof(float));
        uint8_t *out_char = alloca(coded_bits_per_frame);
#else
        COMP  codeword_symbols_de[coded_syms_per_frame];
        float codeword_amps_de[coded_syms_per_frame];
        float llr[coded_bits_per_frame];
        uint8_t out_char[coded_bits_per_frame];
#endif
        gp_deinterleave_comp (codeword_symbols_de, codeword_symbols, coded_syms_per_frame);
        gp_deinterleave_float(codeword_amps_de   , codeword_amps   , coded_syms_per_frame);


        if (f->test_frames) {
            int tmp;
            Nerrs_raw = count_uncoded_errors(ldpc, &f->ofdm->config, &tmp, codeword_symbols_de);
            f->total_bit_errors += Nerrs_raw;
            f->total_bits += f->ofdm_bitsperframe;
        }

        f->modem_frame_count_rx = 0;

        symbols_to_llrs(llr, codeword_symbols_de, codeword_amps_de,
                EsNo, ofdm->mean_amp, coded_syms_per_frame);
        iter = run_ldpc_decoder(ldpc, out_char, llr, &parityCheckCount);

        if (parityCheckCount != ldpc->NumberParityBits)
            rx_status |= RX_BIT_ERRORS;

        if (f->test_frames) {
#ifdef NO_C99
            uint8_t *payload_data_bits = alloca(data_bits_per_frame);
#else
            uint8_t payload_data_bits[data_bits_per_frame];
#endif
            ofdm_generate_payload_data_bits(payload_data_bits, data_bits_per_frame);
            Nerrs_coded = count_errors(payload_data_bits, out_char, data_bits_per_frame);
            f->total_bit_errors_coded += Nerrs_coded;
            f->total_bits_coded += data_bits_per_frame;
        } else {
            memcpy(f->rx_payload_bits, out_char, data_bits_per_frame);
        }

        rx_status |= RX_BITS;

        /* If modem is synced we can decode txt bits */        
        for(k=0; k<f->ofdm_ntxtbits; k++)  { 
            //fprintf(stderr, "txt_bits[%d] = %d\n", k, rx_bits[i]);
            n_ascii = varicode_decode(&f->varicode_dec_states, &ascii_out, &txt_bits[k], 1, 1);
            if (n_ascii && (f->freedv_put_next_rx_char != NULL)) {
                (*f->freedv_put_next_rx_char)(f->callback_state, ascii_out);
            }
        }

        /* estimate uncoded BER from UW */
        for(i=0; i<f->ofdm_nuwbits; i++) {         
            if (rx_uw[i] != ofdm->tx_uw[i]) {
                f->total_bit_errors++;
            }
        }
        
        f->total_bits += f->ofdm_nuwbits;
    }
    
    /* iterate state machine and update nin for next call */
    
    f->nin = ofdm_get_nin(ofdm);
    ofdm_sync_state_machine(ofdm, rx_uw);

    if ((f->verbose && (ofdm->last_sync_state == search)) || (f->verbose == 2)) {
        fprintf(stderr, "%3d nin: %4d st: %-6s euw: %2d %1d f: %5.1f phbw: %d snr: %4.1f eraw: %3d ecdd: %3d iter: %3d pcc: %3d rxst: %s\n",
                f->frames++, ofdm->nin, ofdm_statemode[ofdm->last_sync_state], ofdm->uw_errors, ofdm->sync_counter, 
		(double)ofdm->foff_est_hz, ofdm->phase_est_bandwidth,
                f->snr_est, Nerrs_raw, Nerrs_coded, iter, parityCheckCount, rx_sync_flags_to_text[rx_status]);
    }
    
    return rx_status;
}

