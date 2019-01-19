/*---------------------------------------------------------------------------*\

  FILE........: ofdm_demod.c
  AUTHOR......: David Rowe
  DATE CREATED: Mar 2018

  Given an input file of raw file (8kHz, 16 bit shorts) of OFDM modem
  samples.  Optionally:

    1/ outputs one char per bit (hard decision)
    2/ bit LLRS, one double per bir, for external LDPC decoder like ldpc_dec
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

static struct OFDM_CONFIG *ofdm_config;

static int ofdm_bitsperframe;
static int ofdm_rowsperframe;
static int ofdm_nuwbits;
static int ofdm_ntxtbits;

static const char *progname;

void opt_help() {
    fprintf(stderr, "\nusage: %s [options]\n\n", progname);
    fprintf(stderr, "  Default output file format is one byte per bit hard decision\n\n");
    fprintf(stderr, "  --in          filename  Name of InputModemRawFile\n");
    fprintf(stderr, "  --out         filename  Name of OutputOneCharPerBitFile\n");
    fprintf(stderr, "  --log         filename  Octave log file for testing\n");
    fprintf(stderr, "  --nc          [17..64]  Number of Carriers (17 default, 64 max)\n");
    fprintf(stderr, "  --tcp         Nsecs     Cyclic Prefix Duration (.002 default)\n");
    fprintf(stderr, "  --ts          Nsecs     Symbol Duration (.018 default)\n");
    fprintf(stderr, "  --interleave  depth     Interleaver for LDPC frames, e.g. 1,2,4,8,16 (default is 1)\n");
    fprintf(stderr, "  --tx_freq     freq      Set modulation TX centre Frequency (1500.0 default)\n");
    fprintf(stderr, "  --rx_freq     freq      Set modulation RX centre Frequency (1500.0 default)\n");
    fprintf(stderr, "  --verbose     [1|2]     Verbose output level to stderr (default off)\n");
    fprintf(stderr, "  --testframes            Receive test frames and count errors\n");
    fprintf(stderr, "  --llr                   LLR output boolean, one double per bit\n");
    fprintf(stderr, "  --ldpc                  Run LDPC decoder boolean. This forces 112, one char/bit output values\n");
    fprintf(stderr, "                          per frame.  In testframe mode raw and coded errors will be counted\n\n");

    exit(-1);
}

int main(int argc, char *argv[])
{
    int  i, j, opt;
    char *fin_name, *fout_name, *log_name;

    char *pn = argv[0] + strlen (argv[0]);

    while (pn != argv[0] && !IS_DIR_SEPARATOR (pn[-1]))
        --pn;
    
    progname = pn;

    /* See if they want help */

    if (argc == 1) {
        opt_help();
    }

    /* Turn off stream buffering */

    setvbuf(stdin, NULL, _IONBF, BUFSIZ);
    setvbuf(stdout, NULL, _IONBF, BUFSIZ);

    FILE *fin = stdin;
    FILE *fout = stdout;
    FILE *foct = NULL;

    int input_specified = 0;
    int output_specified = 0;
    int log_specified = 0;
    int logframes = NFRAMES;
    int ldpc_en = 0;
    int llr_en = 0;
    int testframes = 0;
    int interleave_frames = 1;
    int verbose = 0;
    int nc = 17;

    float tcp = 0.002f;
    float ts = 0.018f;
    float rx_centre = 1500.0f;
    float tx_centre = 1500.0f;

    struct optparse options;

    struct optparse_long longopts[] = {
        {"in",         'a', OPTPARSE_REQUIRED},
        {"out",        'b', OPTPARSE_REQUIRED},
        {"log",        'c', OPTPARSE_REQUIRED},
        {"interleave", 'e', OPTPARSE_REQUIRED},
        {"tx_freq",    'f', OPTPARSE_REQUIRED},
        {"rx_freq",    'g', OPTPARSE_REQUIRED},
        {"verbose",    'v', OPTPARSE_REQUIRED},
        {"testframes", 'd', OPTPARSE_NONE},
        {"llr",        'h', OPTPARSE_NONE},
        {"ldpc",       'i', OPTPARSE_NONE},
        {"nc",         'j', OPTPARSE_REQUIRED},
        {"tcp",        'k', OPTPARSE_REQUIRED},
        {"ts",         'l', OPTPARSE_REQUIRED},
        {0, 0, 0}
    };

    optparse_init(&options, argv);

    while ((opt = optparse_long(&options, longopts, NULL)) != -1) {
        switch (opt) {
            case '?':
                opt_help();
            case 'a':
                fin_name = options.optarg;
                input_specified = 1;
                break;
            case 'b':
                fout_name = options.optarg;
                output_specified = 1;
                break;
            case 'c':
                log_name = options.optarg;
                log_specified = 1;
                break;
            case 'd':
                testframes = 1;
                break;
            case 'e':
                interleave_frames = atoi(options.optarg);
            case 'i':                                       /* fall through */
                ldpc_en = 1;
                llr_en = 1;
                break;
            case 'f':
                tx_centre = atof(options.optarg);
                break;
            case 'g':
                rx_centre = atof(options.optarg);
                break;
            case 'h':
                llr_en = 1;
                break;
            case 'j':
                nc = atoi(options.optarg);
                break;
            case 'k':
                tcp = atof(options.optarg);
                break;
            case 'l':
                ts = atof(options.optarg);
                break;
            case 'v':
                verbose = atoi(options.optarg);
                if (verbose < 0 || verbose > 2)
                    verbose = 0;
        }
    }

    /* Print remaining arguments to give user a hint */

    char *arg;

    while ((arg = optparse_arg(&options)))
        fprintf(stderr, "%s\n", arg);

    if (input_specified) {
        if ((fin = fopen(fin_name, "rb")) == NULL) {
            fprintf(stderr, "Error opening input modem sample file: %s\n", fin_name);
            exit(-1);
        }
    }

    if (output_specified) {
        if ((fout = fopen(fout_name, "wb")) == NULL) {
            fprintf(stderr, "Error opening output file: %s\n", fout_name);
            exit(-1);
        }
    }

    if (log_specified) {
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
    ofdm_config->ns = 8;                        /* Number of Symbol frames */
    ofdm_config->bps = 2;                       /* Bits per Symbol */
    ofdm_config->ts = ts;
    ofdm_config->tcp = tcp;
    ofdm_config->tx_centre = tx_centre;
    ofdm_config->rx_centre = rx_centre;
    ofdm_config->fs = 8000.0f;                  /* Sample Frequency */
    ofdm_config->txtbits = 4;                   /* number of auxiliary data bits */
    ofdm_config->state_str = 16;                /* state string length */
    ofdm_config->ftwindowwidth = 11;
    ofdm_config->ofdm_timing_mx_thresh = 0.30f;

    struct OFDM *ofdm = ofdm_create(ofdm_config);
    assert(ofdm != NULL);

    free(ofdm_config);

    /* Get a copy of the actual modem config */
    ofdm_config = ofdm_get_config_param();

    ofdm_bitsperframe = ofdm_get_bits_per_frame();
    ofdm_rowsperframe = ofdm_bitsperframe / (ofdm_config->nc * ofdm_config->bps);
    ofdm_nuwbits = (ofdm_config->ns - 1) * ofdm_config->bps - ofdm_config->txtbits;
    ofdm_ntxtbits = ofdm_config->txtbits;

    float          phase_est_pilot_log[ofdm_rowsperframe * NFRAMES][ofdm_config->nc];
    COMP           rx_np_log[ofdm_rowsperframe * ofdm_config->nc * NFRAMES];
    float          rx_amp_log[ofdm_rowsperframe * ofdm_config->nc * NFRAMES];
    float          foff_hz_log[NFRAMES], snr_est_log[NFRAMES];
    int            timing_est_log[NFRAMES];

    /* zero out the log arrays incase we don't run for NFRAMES and fill them with data */

    for(i=0; i< ofdm_rowsperframe * NFRAMES; i++) {
        for(j=0; j< ofdm_config->nc; j++) {
            phase_est_pilot_log[i][j] = 0.0f;
        }
    }

    for(i=0; i<ofdm_rowsperframe*ofdm_config->nc*NFRAMES; i++) {
        rx_np_log[i].real = 0.0f;
        rx_np_log[i].imag = 0.0f;
        rx_amp_log[i] = 0.0f;
    }

    for(i=0; i<NFRAMES; i++) {
        foff_hz_log[i] = 0.0f;
        snr_est_log[i] = 0.0f;
        timing_est_log[i] = 0.0f;
    }

    /* Set up default LPDC code.  We could add other codes here if we like */

    struct LDPC ldpc;

    set_up_hra_112_112(&ldpc, ofdm_config);

    int data_bits_per_frame = ldpc.data_bits_per_frame;
    int coded_bits_per_frame = ldpc.coded_bits_per_frame;
    int coded_syms_per_frame = ldpc.coded_syms_per_frame;

    if (verbose) {
        fprintf(stderr, "interleave_frames: %d\n",  interleave_frames);
        ofdm_set_verbose(ofdm, verbose);
    }

    int Nerrs_raw[interleave_frames];
    int Nerrs_coded[interleave_frames];
    int iter[interleave_frames];
    int parityCheckCount[interleave_frames];

    for(i=0; i<interleave_frames; i++) {
        Nerrs_raw[i] = 0;
        Nerrs_coded[i] = 0;
        iter[i] = 0;
        parityCheckCount[i] = 0;
    }

    int Nbitsperframe = ofdm_bitsperframe;
    int Nmaxsamperframe = ofdm_get_max_samples_per_frame();

    short  rx_scaled[Nmaxsamperframe];
    int    rx_bits[Nbitsperframe];
    char   rx_bits_char[Nbitsperframe];
    int    rx_uw[ofdm_nuwbits];
    short  txt_bits[ofdm_ntxtbits];
    int    Nerrs, Terrs, Tbits, Terrs2, Tbits2, Terrs_coded, Tbits_coded, frame_count;

    Nerrs = Terrs = Tbits = Terrs2 = Tbits2 = Terrs_coded = Tbits_coded = frame_count = 0;

    float EsNo = 3.0f;

    if (verbose)
        fprintf(stderr,"Warning EsNo: %f hard coded\n", EsNo);

    float snr_est_smoothed_dB = 0.0f;

    COMP  payload_syms[coded_syms_per_frame];
    float payload_amps[coded_syms_per_frame];
    COMP  codeword_symbols[interleave_frames*coded_syms_per_frame];
    float codeword_amps[interleave_frames*coded_syms_per_frame];

    for (i=0; i<interleave_frames*coded_syms_per_frame; i++) {
        codeword_symbols[i].real = 0.0f;
        codeword_symbols[i].imag = 0.0f;
        codeword_amps[i] = 0.0f;
    }

    /* More logging */
    COMP  payload_syms_log[NFRAMES][coded_syms_per_frame];
    float payload_amps_log[NFRAMES][coded_syms_per_frame];

    for (i=0; i<NFRAMES; i++) {
      for (j=0; j<coded_syms_per_frame; j++) {
        payload_syms_log[i][j].real = 0.0f;
        payload_syms_log[i][j].imag = 0.0f;
        payload_amps_log[i][j] = 0.0f;
      }
    }

    int nin_frame = ofdm_get_nin(ofdm);

    int f = 0;

    while(fread(rx_scaled, sizeof(short), nin_frame, fin) == nin_frame) {

        int log_payload_syms = 0;

        /* demod */

        if (strcmp(ofdm->sync_state,"search") == 0) {
            ofdm_sync_search_shorts(ofdm, rx_scaled, (OFDM_AMP_SCALE / 2.0f));
        }

        if ((strcmp(ofdm->sync_state,"synced") == 0) || (strcmp(ofdm->sync_state,"trial") == 0) ) {
            ofdm_demod_shorts(ofdm, rx_bits, rx_scaled, (OFDM_AMP_SCALE / 2.0f));
            ofdm_disassemble_modem_frame(ofdm, rx_uw, payload_syms, payload_amps, txt_bits);
            log_payload_syms = 1;

            /* SNR estimation and smoothing */

            float snr_est_dB = 10.0f * 
                   log10f((ofdm->sig_var/ofdm->noise_var) * ofdm_config->nc * ofdm_config->rs / 3000.0f);

            snr_est_smoothed_dB = 0.9f * snr_est_smoothed_dB + 0.1f * snr_est_dB;

            if (llr_en) {

                /* first few symbols are used for UW and txt bits, find start of (224,112) LDPC codeword
                   and extract QPSK symbols and amplitude estimates */

                assert((ofdm_nuwbits+ofdm_ntxtbits+coded_bits_per_frame) == ofdm_bitsperframe);

                /* now we need to buffer for de-interleaving -------------------------------------*/

                /* shift interleaved symbol buffers to make room for new symbols */

                for(i=0, j=coded_syms_per_frame; j<interleave_frames*coded_syms_per_frame; i++,j++) {
                    codeword_symbols[i] = codeword_symbols[j];
                    codeword_amps[i] = codeword_amps[j];
                }

                /* newest symbols at end of buffer (uses final i from last loop) */

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

                if (ldpc_en) {
                    char out_char[coded_bits_per_frame];

                    interleaver_sync_state_machine(ofdm, &ldpc, ofdm_config, codeword_symbols_de, codeword_amps_de, EsNo,
                                                   interleave_frames, iter, parityCheckCount, Nerrs_coded);

                    if (!strcmp(ofdm->sync_state_interleaver,"synced") && (ofdm->frame_count_interleaver == interleave_frames)) {
                        ofdm->frame_count_interleaver = 0;
                        // if (verbose)
                        //     printf("decode!\n");

                        if (testframes) {
                            Terrs += count_uncoded_errors(&ldpc, ofdm_config, Nerrs_raw, interleave_frames, codeword_symbols_de);
                            Tbits += coded_bits_per_frame*interleave_frames; /* not counting errors in txt bits */
                        }

                        for (j=0; j<interleave_frames; j++) {
                            symbols_to_llrs(llr, &codeword_symbols_de[j*coded_syms_per_frame],
                                                 &codeword_amps_de[j*coded_syms_per_frame],
                                                 EsNo, ofdm->mean_amp, coded_syms_per_frame);
                            iter[j] = run_ldpc_decoder(&ldpc, out_char, llr, &parityCheckCount[j]);

                            //if (verbose)
                            //    fprintf(stderr,"j: %d iter: %d pcc: %d\n", j, iter[j], parityCheckCount[j]);

                            if (testframes) {
                                /* construct payload data bits */
                                
                                int payload_data_bits[data_bits_per_frame];
                                ofdm_generate_payload_data_bits(payload_data_bits, data_bits_per_frame);
                                
                                Nerrs_coded[j] = count_errors(payload_data_bits, out_char, data_bits_per_frame);
                                Terrs_coded += Nerrs_coded[j];
                                Tbits_coded += data_bits_per_frame;
                            }

                            fwrite(out_char, sizeof(char), data_bits_per_frame, fout);
                        }
                    } /* if interleaver synced ..... */

                } else {
                    /* lpdc_en == 0,  external LDPC decoder, so output LLRs */
                    symbols_to_llrs(llr, codeword_symbols_de, codeword_amps_de, EsNo, ofdm->mean_amp, coded_syms_per_frame);
                    fwrite(llr, sizeof(double), coded_bits_per_frame, fout);
                }
            } else {
                /* simple hard decision output for uncoded testing, all bits in frame dumped inlcuding UW and txt */

                for(i=0; i<Nbitsperframe; i++) {
                    rx_bits_char[i] = rx_bits[i];
                }

                fwrite(rx_bits_char, sizeof(char), Nbitsperframe, fout);
            }

            /* optional error counting on uncoded data in non-LDPC testframe mode */

            if (testframes && (ldpc_en == 0)) {
                /* build up a test frame consisting of unique word, txt bits, and psuedo-random
                   uncoded payload bits.  The psuedo-random generator is the same as Octave so
                   it can interoperate with ofdm_tx.m/ofdm_rx.m */

                int Npayloadbits = Nbitsperframe-(ofdm_nuwbits+ofdm_ntxtbits);
                uint16_t r[Npayloadbits];
                uint8_t  payload_bits[Npayloadbits];
                uint8_t  tx_bits[Npayloadbits];
                
                ofdm_rand(r, Npayloadbits);

                for(i=0; i<Npayloadbits; i++) {
                    payload_bits[i] = r[i] > 16384;
                    //if (verbose)
                    //    fprintf(stderr,"%d %d ", r[j], tx_bits_char[i]);
                }

                uint8_t txt_bits[ofdm_ntxtbits];

                for(i=0; i<ofdm_ntxtbits; i++) {
                    txt_bits[i] = 0;
                }

                ofdm_assemble_modem_frame(tx_bits, payload_bits, txt_bits);

                Nerrs = 0;

                for(i=0; i<Nbitsperframe; i++) {
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

        if (ofdm->sync_start) {
            Terrs = Tbits = Terrs2 = Tbits2 = Terrs_coded = Tbits_coded = frame_count = 0;

            for(i=0; i<interleave_frames; i++) {
                Nerrs_raw[i] = Nerrs_coded[i] = 0;
            }

        }

        if (verbose == 2) {
            int r = 0;

            if (testframes) {
                r = (ofdm->frame_count_interleaver - 1 ) % interleave_frames;
            }

            fprintf(stderr, "%3d st: %-6s euw: %2d %1d f: %5.1f ist: %-6s %2d eraw: %3d ecdd: %3d iter: %3d pcc: %3d\n",
                f, ofdm->last_sync_state, ofdm->uw_errors, ofdm->sync_counter, ofdm->foff_est_hz,
                ofdm->last_sync_state_interleaver, ofdm->frame_count_interleaver,
                Nerrs_raw[r], Nerrs_coded[r], iter[r], parityCheckCount[r]);
        }

        /* optional logging of states */

        if (log_specified) {
            /* note corrected phase (rx no phase) is one big linear array for frame */

            for (i = 0; i < ofdm_rowsperframe*ofdm_config->nc; i++) {
                rx_np_log[ofdm_rowsperframe*ofdm_config->nc * f + i].real = crealf(ofdm->rx_np[i]);
                rx_np_log[ofdm_rowsperframe*ofdm_config->nc * f + i].imag = cimagf(ofdm->rx_np[i]);
            }

            /* note phase/amp ests the same for each col, but check them all anyway */

            for (i = 0; i < ofdm_rowsperframe; i++) {
                for (j = 0; j < ofdm_config->nc; j++) {
                    phase_est_pilot_log[ofdm_rowsperframe*f+i][j] = ofdm->aphase_est_pilot_log[ofdm_config->nc*i+j];
                    rx_amp_log[ofdm_rowsperframe*ofdm_config->nc*f+ofdm_config->nc*i+j] = ofdm->rx_amp[ofdm_config->nc*i+j];
                }
            }

            foff_hz_log[f] = ofdm->foff_est_hz;
            timing_est_log[f] = ofdm->timing_est + 1;     /* offset by 1 to match Octave */

            snr_est_log[f] = snr_est_smoothed_dB;

            if (log_payload_syms) {
                for (i=0; i<coded_syms_per_frame; i++) {
                    payload_syms_log[f][i].real = payload_syms[i].real;
                    payload_syms_log[f][i].imag = payload_syms[i].imag;
                    payload_amps_log[f][i] = payload_amps[i];
                }
            }

            if (f == (logframes-1))
                log_specified = 0;
        }

        f++;
    }

    if (input_specified)
        fclose(fin);

    if (output_specified)
        fclose(fout);

    /* optionally dump Octave files */

    if (log_specified) {
        octave_save_float(foct, "phase_est_pilot_log_c", (float*)phase_est_pilot_log, ofdm_rowsperframe*NFRAMES, ofdm_config->nc, ofdm_config->nc);
        octave_save_complex(foct, "rx_np_log_c", (COMP*)rx_np_log, 1, ofdm_rowsperframe*ofdm_config->nc*NFRAMES, ofdm_rowsperframe*ofdm_config->nc*NFRAMES);
        octave_save_float(foct, "rx_amp_log_c", (float*)rx_amp_log, 1, ofdm_rowsperframe*ofdm_config->nc*NFRAMES, ofdm_rowsperframe*ofdm_config->nc*NFRAMES);
        octave_save_float(foct, "foff_hz_log_c", foff_hz_log, NFRAMES, 1, 1);
        octave_save_int(foct, "timing_est_log_c", timing_est_log, NFRAMES, 1);
        octave_save_float(foct, "snr_est_log_c", snr_est_log, NFRAMES, 1, 1);
        octave_save_complex(foct, "payload_syms_log_c", (COMP*)payload_syms_log, NFRAMES, coded_syms_per_frame, coded_syms_per_frame);
        octave_save_float(foct, "payload_amps_log_c", (float*)payload_amps_log, NFRAMES, coded_syms_per_frame, coded_syms_per_frame);

        fclose(foct);
    }

    ofdm_destroy(ofdm);

    if (testframes) {
        float uncoded_ber = (float)Terrs/Tbits;

        if (verbose) {
            fprintf(stderr, "BER......: %5.4f Tbits: %5d Terrs: %5d\n", uncoded_ber, Tbits, Terrs);

            if (!ldpc_en) {
                fprintf(stderr, "BER2.....: %5.4f Tbits: %5d Terrs: %5d\n", (float)Terrs2/Tbits2, Tbits2, Terrs2);
            }
        }

        if (ldpc_en) {
            float coded_ber = (float)Terrs_coded/Tbits_coded;

            if (verbose)
                fprintf(stderr, "Coded BER: %5.4f Tbits: %5d Terrs: %5d\n", coded_ber, Tbits_coded, Terrs_coded);

            /* set return code for Ctest */

            if ((uncoded_ber >= 0.1f) && (coded_ber >= 0.01f))
                return 1;
        }
    }
    
    return 0;
}


