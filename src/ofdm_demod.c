/*---------------------------------------------------------------------------*\

  FILE........: ofdm_demod.c
  AUTHOR......: David Rowe
  DATE CREATED: Mar 2018

  Given an input file of raw file (8kHz, 16 bit shorts) of OFDM modem
  samples.  Optionally:

    1/ outputs one char per bit (hard decision)
    2/ bit LLRS, one double per bit, for external LDPC decoder like ldpc_dec
    3/ LDPC decoded bits, one char per bit

  Also has test frame modes for uncoded and coded operation.

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

#define OPTPARSE_IMPLEMENTATION
#define OPTPARSE_API static
#include "optparse.h"

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <errno.h>

#include "codec2_ofdm.h"
#include "ofdm_internal.h"
#include "octave.h"
#include "mpdecode_core.h"
#include "gp_interleaver.h"
#include "interldpc.h"

#define IS_DIR_SEPARATOR(c)     ((c) == '/')

#define NFRAMES  100               /* just log the first 100 frames          */
#define NDISCARD 20                /* BER2 measure discards first 20 frames   */
#define FS       8000.0f

static struct OFDM_CONFIG *ofdm_config;

static int ofdm_bitsperframe;
static int ofdm_rowsperframe;
static int ofdm_nuwbits;
static int ofdm_ntxtbits;

static const char *progname;

static const char *statemode[] = {
    "search",
    "trial",
    "synced"
};

void opt_help() {
    fprintf(stderr, "\nusage: %s [options]\n\n", progname);
    fprintf(stderr, "  Default output file format is one byte per bit hard decision\n\n");
    fprintf(stderr, "  --in          filename   Name of InputModemRawFile\n");
    fprintf(stderr, "  --out         filename   Name of OutputOneCharPerBitFile\n");
    fprintf(stderr, "  --log         filename   Octave log file for testing\n");
    fprintf(stderr, "  --nc          [17..62]   Number of Carriers (17 default, 62 max)\n");
    fprintf(stderr, "  --ns           Nframes   Number of Symbol Frames (8 default)\n");
    fprintf(stderr, "  --tcp            Nsecs   Cyclic Prefix Duration (.002 default)\n");
    fprintf(stderr, "  --ts             Nsecs   Symbol Duration (.018 default)\n");
    fprintf(stderr, "  --bandwidth      [0|1]   Select phase est bw mode AUTO low or high (0) or LOCKED high (1) (default 0)\n");
    fprintf(stderr, "  --interleave     depth   Interleaver for LDPC frames, e.g. 1,2,4,8,16 (default is 1)\n");
    fprintf(stderr, "                           Must also specify --ldpc option\n");
    fprintf(stderr, "  --tx_freq         freq   Set modulation TX centre Frequency (1500.0 default)\n");
    fprintf(stderr, "  --rx_freq         freq   Set modulation RX centre Frequency (1500.0 default)\n");
    fprintf(stderr, "  --verbose      [1|2|3]   Verbose output level to stderr (default off)\n");
    fprintf(stderr, "  --testframes             Receive test frames and count errors\n");
    fprintf(stderr, "  --llr                    LLR output boolean, one double per bit\n");
    fprintf(stderr, "  --ldpc           [1|2]   Run LDPC decoder In (224,112) 700D or (504, 396) 2020 mode.\n");
    fprintf(stderr, "  --databits     numBits   Number of data bits used in LDPC codeword.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  --start_secs      secs   Number of seconds delay before we start to demod\n");
    fprintf(stderr, "  --len_secs        secs   Number of seconds to run demod\n");
    fprintf(stderr, "  --skip_secs   timeSecs   At timeSecs introduce a large timing error by skipping half a frame of samples\n");
    fprintf(stderr, "  --dpsk                   Differential PSK.\n");
    fprintf(stderr, "\n");
    
    exit(-1);
}

int main(int argc, char *argv[]) {
    int i, j, opt, val;

    char *pn = argv[0] + strlen(argv[0]);

    while (pn != argv[0] && !IS_DIR_SEPARATOR(pn[-1]))
        --pn;

    progname = pn;

    /* Turn off stream buffering */

    setvbuf(stdin, NULL, _IONBF, BUFSIZ);
    setvbuf(stdout, NULL, _IONBF, BUFSIZ);

    FILE *fin = stdin;
    FILE *fout = stdout;
    FILE *foct = NULL;

    char *fin_name = NULL;
    char *fout_name = NULL;
    char *log_name = NULL;

    int logframes = NFRAMES;
    int interleave_frames = 1;
    int nc = 17;
    int ns = 8;
    int verbose = 0;
    int phase_est_bandwidth_mode = AUTO_PHASE_EST;
    int ldpc_en = 0;
    int data_bits_per_frame = 0;

    bool testframes = false;
    bool input_specified = false;
    bool output_specified = false;
    bool log_specified = false;
    bool log_active = false;
    bool ldpc_check = false;
    bool llr_en = false;
    bool dpsk = false;
    
    float tcp = 0.002f;
    float ts = 0.018f;
    float rx_centre = 1500.0f;
    float tx_centre = 1500.0f;

    float time_to_sync = -1;
    float start_secs = 0.0;
    float len_secs = 0.0;
    float skip_secs = 0.0;
    struct optparse options;

    struct optparse_long longopts[] = {
        {"in", 'a', OPTPARSE_REQUIRED},
        {"out", 'b', OPTPARSE_REQUIRED},
        {"log", 'c', OPTPARSE_REQUIRED},
        {"testframes", 'd', OPTPARSE_NONE},
        {"interleave", 'e', OPTPARSE_REQUIRED},
        {"bandwidth", 'o', OPTPARSE_REQUIRED},
        {"tx_freq", 'f', OPTPARSE_REQUIRED},
        {"rx_freq", 'g', OPTPARSE_REQUIRED},
        {"verbose", 'v', OPTPARSE_REQUIRED},
        {"llr", 'h', OPTPARSE_NONE},
        {"ldpc", 'i', OPTPARSE_REQUIRED},
        {"nc", 'j', OPTPARSE_REQUIRED},
        {"tcp", 'k', OPTPARSE_REQUIRED},
        {"ts", 'l', OPTPARSE_REQUIRED},
        {"ns", 'm', OPTPARSE_REQUIRED},
        {"databits", 'p', OPTPARSE_REQUIRED},        
        {"start_secs", 'x', OPTPARSE_REQUIRED},        
        {"len_secs", 'y', OPTPARSE_REQUIRED},        
        {"skip_secs", 'z', OPTPARSE_REQUIRED},        
        {"dpsk", 'q', OPTPARSE_NONE},        
        {0, 0, 0}
    };

    optparse_init(&options, argv);

    while ((opt = optparse_long(&options, longopts, NULL)) != -1) {
        switch (opt) {
            case '?':
                opt_help();
            case 'a':
                fin_name = options.optarg;
                input_specified = true;
                break;
            case 'b':
                fout_name = options.optarg;
                output_specified = true;
                break;
            case 'c':
                log_name = options.optarg;
                log_specified = true;
                log_active = true;
                break;
            case 'd':
                testframes = true;
                break;
            case 'e':
                interleave_frames = atoi(options.optarg);
                ldpc_check = true;
                llr_en = true;
                break;
            case 'i':
                ldpc_en = atoi(options.optarg);
                if ((ldpc_en != 1) && (ldpc_en !=2)) {
                    fprintf(stderr, "--ldpc 1  (224,112) code used for 700D\n");
                    fprintf(stderr, "--ldpc 2  (504,396) code used for 2020\n");
                    opt_help();
                }
            case 'h': /* fall through */
                llr_en = true;
                break;
            case 'f':
                tx_centre = atof(options.optarg);
                break;
            case 'g':
                rx_centre = atof(options.optarg);
                break;
            case 'j':
                val = atoi(options.optarg);

                if (val > 62 || val < 17) {
                    opt_help();
                } else {
                    nc = val;
                }
                break;
            case 'k':
                tcp = atof(options.optarg);
                break;
            case 'l':
                ts = atof(options.optarg);
                break;
            case 'm':
                ns = atoi(options.optarg);
                break;
            case 'o':
                phase_est_bandwidth_mode = atoi(options.optarg);
                break;
            case 'p':
                data_bits_per_frame = atoi(options.optarg);
                break;
            case 'q':
                dpsk = true;
                break;
            case 'v':
                verbose = atoi(options.optarg);
                if (verbose < 0 || verbose > 3)
                    verbose = 0;
                break;
            case 'x':
                 start_secs = atoi(options.optarg);
                 break;
            case 'y':
                 len_secs = atoi(options.optarg);
                 break;
            case 'z':
                 skip_secs = atoi(options.optarg);
                 break;
                 
       }
    }

    /* Print remaining arguments to give user a hint */
    char *arg;

    while ((arg = optparse_arg(&options)))
        fprintf(stderr, "%s\n", arg);

    if (ldpc_check == true) {
        if ((interleave_frames != 1) && (ldpc_en == 0)) {
            fprintf(stderr, "You must also specify an ldpc option with interleave\n");
            exit(-1);
        }
    }

    if (input_specified == true) {
        if ((fin = fopen(fin_name, "rb")) == NULL) {
            fprintf(stderr, "Error opening input modem sample file: %s\n", fin_name);
            exit(-1);
        }
    }

    if (output_specified == true) {
        if ((fout = fopen(fout_name, "wb")) == NULL) {
            fprintf(stderr, "Error opening output file: %s\n", fout_name);
            exit(-1);
        }
    }

    if (log_specified == true) {
        if ((foct = fopen(log_name, "wt")) == NULL) {
            fprintf(stderr, "Error opening Octave output file: %s\n", log_name);
            exit(-1);
        }
    }

    if ((ofdm_config = (struct OFDM_CONFIG *) calloc(1, sizeof (struct OFDM_CONFIG))) == NULL) {
        fprintf(stderr, "Out of Memory\n");
        exit(-1);
    }

    ofdm_config->nc = nc;
    ofdm_config->ns = ns; /* Number of Symbol frames */
    ofdm_config->bps = 2; /* Bits per Symbol */
    ofdm_config->ts = ts;
    ofdm_config->tcp = tcp;
    ofdm_config->tx_centre = tx_centre;
    ofdm_config->rx_centre = rx_centre;
    ofdm_config->fs = FS; /* Sample Frequency */
    ofdm_config->txtbits = 4; /* number of auxiliary data bits */

    ofdm_config->ftwindowwidth = 11;
    ofdm_config->ofdm_timing_mx_thresh = 0.30f;

    struct OFDM *ofdm = ofdm_create(ofdm_config);
    assert(ofdm != NULL);

    free(ofdm_config);

    ofdm_set_phase_est_bandwidth_mode(ofdm, phase_est_bandwidth_mode);
    ofdm_set_dpsk(ofdm, dpsk);
    
    /* Get a copy of the actual modem config */
    ofdm_config = ofdm_get_config_param();

    ofdm_bitsperframe = ofdm_get_bits_per_frame();
    ofdm_rowsperframe = ofdm_bitsperframe / (ofdm_config->nc * ofdm_config->bps);
    ofdm_nuwbits = (ofdm_config->ns - 1) * ofdm_config->bps - ofdm_config->txtbits;
    ofdm_ntxtbits = ofdm_config->txtbits;

    float phase_est_pilot_log[ofdm_rowsperframe * NFRAMES][ofdm_config->nc];
    COMP rx_np_log[ofdm_rowsperframe * ofdm_config->nc * NFRAMES];
    float rx_amp_log[ofdm_rowsperframe * ofdm_config->nc * NFRAMES];
    float foff_hz_log[NFRAMES], snr_est_log[NFRAMES];
    int timing_est_log[NFRAMES];

    /* zero out the log arrays in case we don't run for NFRAMES and fill them with data */

    for (i = 0; i < (ofdm_rowsperframe * NFRAMES); i++) {
        for (j = 0; j < ofdm_config->nc; j++) {
            phase_est_pilot_log[i][j] = 0.0f;
        }
    }

    for (i = 0; i < (ofdm_rowsperframe * ofdm_config->nc * NFRAMES); i++) {
        rx_np_log[i].real = 0.0f;
        rx_np_log[i].imag = 0.0f;
        rx_amp_log[i] = 0.0f;
    }

    for (i = 0; i < NFRAMES; i++) {
        foff_hz_log[i] = 0.0f;
        snr_est_log[i] = 0.0f;
        timing_est_log[i] = 0.0f;
    }

    /* Set up default LPDC code.  We could add other codes here if we like */

    struct LDPC ldpc;

    int coded_bits_per_frame = 0;
    int coded_syms_per_frame = 0;
    
    if (ldpc_en) {
        if (ldpc_en == 1)
            set_up_hra_112_112(&ldpc, ofdm_config);
        else
            set_up_hra_504_396(&ldpc, ofdm_config);

        /* here is where we can change data bits per frame to a number smaller than LDPC code input data bits_per_frame */
        if (data_bits_per_frame) {
            set_data_bits_per_frame(&ldpc, data_bits_per_frame, ofdm_config->bps);
        }
    
        data_bits_per_frame = ldpc.data_bits_per_frame;
        coded_bits_per_frame = ldpc.coded_bits_per_frame;
        coded_syms_per_frame = ldpc.coded_syms_per_frame;
 
        assert(data_bits_per_frame <= ldpc.ldpc_data_bits_per_frame);
        assert(coded_bits_per_frame <= ldpc.ldpc_coded_bits_per_frame);
        
        if (verbose > 1) {
            fprintf(stderr, "ldpc_data_bits_per_frame = %d\n", ldpc.ldpc_data_bits_per_frame);
            fprintf(stderr, "ldpc_coded_bits_per_frame  = %d\n", ldpc.ldpc_coded_bits_per_frame);
            fprintf(stderr, "data_bits_per_frame = %d\n", data_bits_per_frame);
            fprintf(stderr, "coded_bits_per_frame  = %d\n", coded_bits_per_frame);
            fprintf(stderr, "ofdm_bits_per_frame  = %d\n", ofdm_bitsperframe);
        }

    }
    
    if (verbose != 0) {
        fprintf(stderr, "interleave_frames: %d\n", interleave_frames);
        ofdm_set_verbose(ofdm, verbose);
        
        fprintf(stderr, "Phase Estimate Switching: ");

        switch (phase_est_bandwidth_mode) {
        case 0: fprintf(stderr, "Auto\n");
                break;
        case 1: fprintf(stderr, "Locked\n");
        }
    }

    int Nerrs_raw[interleave_frames];
    int Nerrs_coded[interleave_frames];
    int iter[interleave_frames];
    int parityCheckCount[interleave_frames];

    for (i = 0; i < interleave_frames; i++) {
        Nerrs_raw[i] = 0;
        Nerrs_coded[i] = 0;
        iter[i] = 0;
        parityCheckCount[i] = 0;
    }

    int Nbitsperframe = ofdm_bitsperframe;
    int Nmaxsamperframe = ofdm_get_max_samples_per_frame();
    // TODO: these constants come up a lot so might be best placed in ofdm_create()
    int Npayloadbitsperframe = ofdm_bitsperframe - ofdm_nuwbits - ofdm_ntxtbits;
    int Npayloadsymsperframe = Npayloadbitsperframe/ofdm_config->bps;

    if (ldpc_en) assert(Npayloadsymsperframe >= coded_syms_per_frame);

    short rx_scaled[Nmaxsamperframe];
    int rx_bits[Nbitsperframe];
    uint8_t rx_bits_char[Nbitsperframe];
    uint8_t rx_uw[ofdm_nuwbits];
    short txt_bits[ofdm_ntxtbits];
    int Nerrs, Terrs, Tbits, Terrs2, Tbits2, Terrs_coded, Tbits_coded, frame_count;

    Nerrs = Terrs = Tbits = Terrs2 = Tbits2 = Terrs_coded = Tbits_coded = frame_count = 0;

    float EsNo = 3.0f;
    float snr_est_smoothed_dB = 0.0f;
    
    if (verbose == 2)
        fprintf(stderr, "Warning EsNo: %f hard coded\n", EsNo);

    COMP payload_syms[Npayloadsymsperframe];
    COMP codeword_symbols[interleave_frames * Npayloadsymsperframe];

    float payload_amps[Npayloadsymsperframe];
    float codeword_amps[interleave_frames * Npayloadsymsperframe];

    for (i = 0; i < (interleave_frames * Npayloadsymsperframe); i++) {
        codeword_symbols[i].real = 0.0f;
        codeword_symbols[i].imag = 0.0f;
        codeword_amps[i] = 0.0f;
    }

    /* More logging */
    COMP payload_syms_log[NFRAMES][Npayloadsymsperframe];
    float payload_amps_log[NFRAMES][Npayloadsymsperframe];

    for (i = 0; i < NFRAMES; i++) {
        for (j = 0; j < Npayloadsymsperframe; j++) {
            payload_syms_log[i][j].real = 0.0f;
            payload_syms_log[i][j].imag = 0.0f;
            payload_amps_log[i][j] = 0.0f;
        }
    }

    int nin_frame = ofdm_get_nin(ofdm);

    int f = 0;
    int finish = 0;

    if (start_secs != 0.0) {
        int offset = start_secs*FS*sizeof(short);
        fseek(fin, offset, SEEK_SET);
    }
    
    while ((fread(rx_scaled, sizeof (short), nin_frame, fin) == nin_frame) && !finish) {

        bool log_payload_syms = false;

        /* demod */

        if (ofdm->sync_state == search) {
            ofdm_sync_search_shorts(ofdm, rx_scaled, (OFDM_AMP_SCALE / 2.0f));
        }

        if ((ofdm->sync_state == synced) || (ofdm->sync_state == trial)) {
            ofdm_demod_shorts(ofdm, rx_bits, rx_scaled, (OFDM_AMP_SCALE / 2.0f));
            ofdm_disassemble_modem_frame(ofdm, rx_uw, payload_syms, payload_amps, txt_bits);
            log_payload_syms = true;

            /* SNR estimation and smoothing */

            float snr_est_dB = 10.0f *
                    log10f((ofdm->sig_var / ofdm->noise_var) *
                    ofdm_config->nc * ofdm_config->rs / 3000.0f);

            snr_est_smoothed_dB = 0.9f * snr_est_smoothed_dB + 0.1f * snr_est_dB;

            if (llr_en == true) {

                /* first few symbols are used for UW and txt bits, find start of (224,112) LDPC codeword
                   and extract QPSK symbols and amplitude estimates */

                assert((ofdm_nuwbits + ofdm_ntxtbits + coded_bits_per_frame) <= ofdm_bitsperframe);

                /* now we need to buffer for de-interleaving -------------------------------------*/

                /* shift interleaved symbol buffers to make room for new symbols */

                for (i = 0, j = coded_syms_per_frame; j < interleave_frames * coded_syms_per_frame; i++, j++) {
                    codeword_symbols[i] = codeword_symbols[j];
                    codeword_amps[i] = codeword_amps[j];
                }

                /* newest symbols at end of buffer (uses final i from last loop) */

                for (i = (interleave_frames - 1) * coded_syms_per_frame, j = 0; i < interleave_frames * coded_syms_per_frame; i++, j++) {
                    codeword_symbols[i] = payload_syms[j];
                    codeword_amps[i] = payload_amps[j];
                }

                /* run de-interleaver */

                COMP codeword_symbols_de[interleave_frames * coded_syms_per_frame];
                float codeword_amps_de[interleave_frames * coded_syms_per_frame];

                gp_deinterleave_comp(codeword_symbols_de, codeword_symbols, interleave_frames * coded_syms_per_frame);
                gp_deinterleave_float(codeword_amps_de, codeword_amps, interleave_frames * coded_syms_per_frame);

                float llr[coded_bits_per_frame];

                if (ldpc_en) {
                    uint8_t out_char[coded_bits_per_frame];

                    interleaver_sync_state_machine(ofdm, &ldpc, ofdm_config, codeword_symbols_de, codeword_amps_de, EsNo,
                            interleave_frames, iter, parityCheckCount, Nerrs_coded);

                    if ((ofdm->sync_state_interleaver == synced) && (ofdm->frame_count_interleaver == interleave_frames)) {
                        ofdm->frame_count_interleaver = 0;

                        if (testframes == true) {
                            Terrs += count_uncoded_errors(&ldpc, ofdm_config, Nerrs_raw, interleave_frames, codeword_symbols_de);
                            Tbits += coded_bits_per_frame*interleave_frames; /* not counting errors in txt bits */
                        }

                        for (j = 0; j < interleave_frames; j++) {
                            symbols_to_llrs(llr, &codeword_symbols_de[j * coded_syms_per_frame],
                                    &codeword_amps_de[j * coded_syms_per_frame],
                                    EsNo, ofdm->mean_amp, coded_syms_per_frame);                           
                            if (ldpc.data_bits_per_frame == ldpc.ldpc_data_bits_per_frame) {
                                /* all data bits in code word used */
                               iter[j] = run_ldpc_decoder(&ldpc, out_char, llr, &parityCheckCount[j]);
                            } else {
                                /* some unused data bits, set these to known values to strengthen code */
                                float llr_full_codeword[ldpc.ldpc_coded_bits_per_frame];
                                int unused_data_bits = ldpc.ldpc_data_bits_per_frame - ldpc.data_bits_per_frame;
                                
                                // received data bits
                                for (i = 0; i < ldpc.data_bits_per_frame; i++)
                                    llr_full_codeword[i] = llr[i]; 
                                // known bits ... so really likely
                                for (i = ldpc.data_bits_per_frame; i < ldpc.ldpc_data_bits_per_frame; i++)
                                    llr_full_codeword[i] = -100.0;
                                // parity bits at end
                                for (i = ldpc.ldpc_data_bits_per_frame; i < ldpc.ldpc_coded_bits_per_frame; i++)
                                    llr_full_codeword[i] = llr[i-unused_data_bits]; 
                                iter[j] = run_ldpc_decoder(&ldpc, out_char, llr_full_codeword, &parityCheckCount[j]);
                            }

                            if (testframes == true) {
                                /* construct payload data bits */

                                uint8_t payload_data_bits[data_bits_per_frame];
                                ofdm_generate_payload_data_bits(payload_data_bits, data_bits_per_frame);

                                Nerrs_coded[j] = count_errors(payload_data_bits, out_char, data_bits_per_frame);
                                Terrs_coded += Nerrs_coded[j];
                                Tbits_coded += data_bits_per_frame;
                            }

                            fwrite(out_char, sizeof (char), data_bits_per_frame, fout);
                        }
                    } /* if interleaver synced ..... */
                } else {
                    /* lpdc_en == false,  external LDPC decoder, so output LLRs */
                    symbols_to_llrs(llr, codeword_symbols_de, codeword_amps_de, EsNo, ofdm->mean_amp, coded_syms_per_frame);
                    fwrite(llr, sizeof (double), coded_bits_per_frame, fout);
                }
            } else {
                /* simple hard decision output for uncoded testing, all bits in frame dumped inlcuding UW and txt */

                for (i = 0; i < Nbitsperframe; i++) {
                    rx_bits_char[i] = rx_bits[i];
                }

                fwrite(rx_bits_char, sizeof (uint8_t), Nbitsperframe, fout);
            }

            /* optional error counting on uncoded data in non-LDPC testframe mode */

            if ((testframes == true) && (ldpc_en ==0)) {
                /* build up a test frame consisting of unique word, txt bits, and psuedo-random
                   uncoded payload bits.  The psuedo-random generator is the same as Octave so
                   it can interoperate with ofdm_tx.m/ofdm_rx.m */

                int Npayloadbits = Nbitsperframe - (ofdm_nuwbits + ofdm_ntxtbits);
                uint16_t r[Npayloadbits];
                uint8_t payload_bits[Npayloadbits];
                uint8_t tx_bits[Npayloadbits];

                ofdm_rand(r, Npayloadbits);

                for (i = 0; i < Npayloadbits; i++) {
                    payload_bits[i] = r[i] > 16384;
                }

                uint8_t txt_bits[ofdm_ntxtbits];

                for (i = 0; i < ofdm_ntxtbits; i++) {
                    txt_bits[i] = 0;
                }

                ofdm_assemble_modem_frame(ofdm, tx_bits, payload_bits, txt_bits);

                Nerrs = 0;

                for (i = 0; i < Nbitsperframe; i++) {
                    if (tx_bits[i] != rx_bits[i]) {
                        Nerrs++;
                    }
                }

                Terrs += Nerrs;
                Tbits += Nbitsperframe;

                if (frame_count >= NDISCARD) {
                    Terrs2 += Nerrs;
                    Tbits2 += Nbitsperframe;
                }
            }

            frame_count++;
        }

        nin_frame = ofdm_get_nin(ofdm);
        ofdm_sync_state_machine(ofdm, rx_uw);

        /* act on any events returned by state machine */

        if (ofdm->sync_start == true) {
            Terrs = Tbits = Terrs2 = Tbits2 = Terrs_coded = Tbits_coded = frame_count = 0;

            for (i = 0; i < interleave_frames; i++) {
                Nerrs_raw[i] = 0;
                Nerrs_coded[i] = 0;
            }
        }

        #ifdef TT
        if (ldpc_en) {
            /* experimental timing gear switching based on LDPC TODO: move into shared function */
            bool ldpc_decode_ok = parityCheckCount[0] > 0.8*ldpc.NumberParityBits;
            if (ofdm->sync_state == synced) {
                if ((ofdm->sync_counter>2) && (ldpc_decode_ok == false)) {
                    ofdm_set_timing_range(ofdm, WIDE_TIMING);
                }
                else {
                    ofdm_set_timing_range(ofdm, NARROW_TIMING);
                }
            }
        }
        #endif
        
        if (verbose >= 2) {
            int r = 0;

            if (testframes == true) {
                r = (ofdm->frame_count_interleaver - 1) % interleave_frames;
            }

           fprintf(stderr, "%3d nin: %4d st: %-6s euw: %2d %1d f: %5.1f pbw: %d ist: %-6s %2d eraw: %3d ecdd: %3d iter: %3d pcc: %3d\n",
                    f, nin_frame, 
                    statemode[ofdm->last_sync_state],
                    ofdm->uw_errors,
                    ofdm->sync_counter,
                    ofdm->foff_est_hz,
                    ofdm->phase_est_bandwidth,
                    statemode[ofdm->last_sync_state_interleaver],
                    ofdm->frame_count_interleaver,
                    Nerrs_raw[r], Nerrs_coded[r], iter[r], parityCheckCount[r]);

            /* detect a sucessful sync for time to sync tests */
            if ((time_to_sync < 0) && ((ofdm->sync_state == synced) || (ofdm->sync_state == trial)))          
                if ((parityCheckCount[r] > 80) && (iter[r] != 100))
                    time_to_sync = (float)(f+1)*ofdm_get_samples_per_frame()/FS;

        }

        /* optional logging of states */

        if (log_active == true) {
            /* note corrected phase (rx no phase) is one big linear array for frame */

            for (i = 0; i < ofdm_rowsperframe * ofdm_config->nc; i++) {
                rx_np_log[ofdm_rowsperframe * ofdm_config->nc * f + i].real = crealf(ofdm->rx_np[i]);
                rx_np_log[ofdm_rowsperframe * ofdm_config->nc * f + i].imag = cimagf(ofdm->rx_np[i]);
            }

            /* note phase/amp ests the same for each col, but check them all anyway */

            for (i = 0; i < ofdm_rowsperframe; i++) {
                for (j = 0; j < ofdm_config->nc; j++) {
                    phase_est_pilot_log[ofdm_rowsperframe * f + i][j] = ofdm->aphase_est_pilot_log[ofdm_config->nc * i + j];
                    rx_amp_log[ofdm_rowsperframe * ofdm_config->nc * f + ofdm_config->nc * i + j] = ofdm->rx_amp[ofdm_config->nc * i + j];
                }
            }

            foff_hz_log[f] = ofdm->foff_est_hz;
            timing_est_log[f] = ofdm->timing_est + 1; /* offset by 1 to match Octave */

            snr_est_log[f] = snr_est_smoothed_dB;

            if (log_payload_syms == true) {
                for (i = 0; i < coded_syms_per_frame; i++) {
                    payload_syms_log[f][i].real = payload_syms[i].real;
                    payload_syms_log[f][i].imag = payload_syms[i].imag;
                    payload_amps_log[f][i] = payload_amps[i];
                }
            }

            if (f == (logframes - 1))
                log_active = false;
        }

        if (len_secs != 0.0) {
            float secs = (float)f*ofdm_get_samples_per_frame()/FS;
            if (secs >= len_secs) finish = 1;
        }
        
        if (skip_secs != 0.0) {
            /* big nasty timing error */
            float secs = (float)f*ofdm_get_samples_per_frame()/FS;
            if (secs >= skip_secs) {
                assert(fread(rx_scaled, sizeof (short), nin_frame/2, fin) == nin_frame/2);
                fprintf(stderr,"  Skip!  Just introduced a nasty big timing slip\n");
                skip_secs = 0.0; /* make sure we just introduce one error */
            }
        }
         
        f++;
    }

    ofdm_destroy(ofdm);

    if (input_specified == true)
        fclose(fin);

    if (output_specified == true)
        fclose(fout);

    /* optionally dump Octave files */

    if (log_specified == true) {
        octave_save_float(foct, "phase_est_pilot_log_c", (float*) phase_est_pilot_log, ofdm_rowsperframe*NFRAMES, ofdm_config->nc, ofdm_config->nc);
        octave_save_complex(foct, "rx_np_log_c", (COMP*) rx_np_log, 1, ofdm_rowsperframe * ofdm_config->nc*NFRAMES, ofdm_rowsperframe * ofdm_config->nc * NFRAMES);
        octave_save_float(foct, "rx_amp_log_c", (float*) rx_amp_log, 1, ofdm_rowsperframe * ofdm_config->nc*NFRAMES, ofdm_rowsperframe * ofdm_config->nc * NFRAMES);
        octave_save_float(foct, "foff_hz_log_c", foff_hz_log, NFRAMES, 1, 1);
        octave_save_int(foct, "timing_est_log_c", timing_est_log, NFRAMES, 1);
        octave_save_float(foct, "snr_est_log_c", snr_est_log, NFRAMES, 1, 1);
        octave_save_complex(foct, "payload_syms_log_c", (COMP*) payload_syms_log, NFRAMES, coded_syms_per_frame, coded_syms_per_frame);
        octave_save_float(foct, "payload_amps_log_c", (float*) payload_amps_log, NFRAMES, coded_syms_per_frame, coded_syms_per_frame);

        fclose(foct);
    }

    if (verbose == 2)
        printf("time_to_sync: %f\n", time_to_sync);
    
    if (testframes == true) {
        float uncoded_ber = (float) Terrs / Tbits;

        if (verbose != 0) {
            fprintf(stderr, "BER......: %5.4f Tbits: %5d Terrs: %5d\n", uncoded_ber, Tbits, Terrs);

            if ((ldpc_en == 0) && (frame_count > NDISCARD)) {
                fprintf(stderr, "BER2.....: %5.4f Tbits: %5d Terrs: %5d\n", (float) Terrs2 / Tbits2, Tbits2, Terrs2);
            }
        }

        if (ldpc_en) {
            float coded_ber = (float) Terrs_coded / Tbits_coded;

            if (verbose != 0)
                fprintf(stderr, "Coded BER: %5.4f Tbits: %5d Terrs: %5d\n", coded_ber, Tbits_coded, Terrs_coded);

            /* set return code for Ctest, 1 for fail */

            if ((Tbits == 0) || (Tbits_coded == 0) || (uncoded_ber >= 0.1f) || (coded_ber >= 0.01f))
                return 1;
        }
    }

    return 0;
}
