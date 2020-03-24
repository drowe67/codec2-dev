/*---------------------------------------------------------------------------*\

  FILE........: ofdm_mod.c
  AUTHOR......: David Rowe
  DATE CREATED: March 2018

  Given an input file of bits (note one bit per char format), outputs
  a raw file (8kHz, 16 bit shorts) of OFDM modem samples ready to send
  over a HF radio channel.

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

#include "codec2_ofdm.h"
#include "ofdm_internal.h"
#include "interldpc.h"
#include "gp_interleaver.h"
#include "varicode.h"

#define IS_DIR_SEPARATOR(c)     ((c) == '/')

static struct OFDM_CONFIG *ofdm_config;

static int ofdm_bitsperframe;
static int ofdm_nuwbits;
static int ofdm_ntxtbits;

static const char *progname;

void opt_help() {
    fprintf(stderr, "\nusage: %s [options]\n\n", progname);
    fprintf(stderr, "  --in      filename    Name of InputOneCharPerBitFile\n");
    fprintf(stderr, "  --out     filename    Name of OutputModemRawFile\n");
    fprintf(stderr, "  --nc      [17..62]    Number of Carriers (17 default, 62 max)\n");
    fprintf(stderr, "  --ns       Nframes    Number of Symbol Frames (8 default)\n");
    fprintf(stderr, "  --tcp        Nsecs    Cyclic Prefix Duration (.002 default)\n");
    fprintf(stderr, "  --ts         Nsecs    Symbol Duration (.018 default)\n");
    fprintf(stderr, "  --testframes Nsecs    Transmit test frames (adjusts test frames for raw and LDPC modes)\n");
    fprintf(stderr, "  --interleave depth    Interleave depth for LDPC frames, e.g. 1,2,4,8,16 (default is 1)\n");
    fprintf(stderr, "  --tx_freq     freq    Set an optional modulation TX centre frequency (1500.0 default)\n");
    fprintf(stderr, "  --rx_freq     freq    Set an optional modulation RX centre frequency (1500.0 default)\n\n");
    fprintf(stderr, "  --verbose  [1|2|3]    Verbose output level to stderr (default off)\n");
    fprintf(stderr, "  --txbpf               Transmit band pass filter boolean (default off)\n");
    fprintf(stderr, "  --text                Include a standard text message boolean (default off)\n");
    fprintf(stderr, "  -i --ldpc    [1|2]    Run LDPC decoder (1 -> (224,112) 700D code, 2 -> (504,396) 2020 code).\n"
                    "                        In testframe mode raw and coded errors will be counted.\n");
    fprintf(stderr, "  -p --databits numBits Number of data bits used in LDPC codeword.\n");
    fprintf(stderr, "  --dpsk                Differential PSK.\n");
    fprintf(stderr, "\n");
    exit(-1);
}

int main(int argc, char *argv[]) {
    char *fin_name, *fout_name;
    int i, j, opt, val;

    char *pn = argv[0] + strlen(argv[0]);

    while (pn != argv[0] && !IS_DIR_SEPARATOR(pn[-1]))
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

    /* set for LDPC coded or uncoded frames */

    int ldpc_en = 0;
    int interleave_frames = 1;

    int input_specified = 0;
    int output_specified = 0;
    int verbose = 0;
    int txbpf_en = 0;
    int testframes = 0;
    int use_text = 0;
    int dpsk = 0;

    int Nframes = 0;
    int Nsec = 0;
    int Nrows = 0;

    int nc = 17;
    int ns = 8;
    float tcp = 0.0020f;
    float ts = 0.0180f;
    float rx_centre = 1500.0f;
    float tx_centre = 1500.0f;

    int   data_bits_per_frame = 0;
    struct optparse options;

    struct optparse_long longopts[] = {
        {"in", 'a', OPTPARSE_REQUIRED},
        {"out", 'b', OPTPARSE_REQUIRED},
        {"nc", 'c', OPTPARSE_REQUIRED},
        {"ns", 'm', OPTPARSE_REQUIRED},
        {"tcp", 'd', OPTPARSE_REQUIRED},
        {"ts", 'e', OPTPARSE_REQUIRED},
        {"testframes", 'f', OPTPARSE_REQUIRED},
        {"interleave", 'g', OPTPARSE_REQUIRED},
        {"tx_freq", 'h', OPTPARSE_REQUIRED},
        {"rx_freq", 'i', OPTPARSE_REQUIRED},
        {"ldpc", 'j', OPTPARSE_REQUIRED},
        {"txbpf", 'k', OPTPARSE_NONE},
        {"text", 'l', OPTPARSE_NONE},
        {"verbose", 'v', OPTPARSE_REQUIRED},
        {"databits", 'p', OPTPARSE_REQUIRED},        
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
                input_specified = 1;
                break;
            case 'b':
                fout_name = options.optarg;
                output_specified = 1;
                break;
            case 'c':
                val = atoi(options.optarg);

                if (val > 62 || val < 17) {
                    opt_help();
                } else {
                    nc = val;
                }
                break;
            case 'd':
                tcp = atof(options.optarg);
                break;
            case 'e':
                ts = atof(options.optarg);
                break;
            case 'm':
                ns = atoi(options.optarg);
                break;
            case 'f':
                testframes = 1;
                Nsec = atoi(options.optarg);
                break;
            case 'g':
                interleave_frames = atoi(options.optarg);
                break;
            case 'h':
                tx_centre = atof(options.optarg);
                break;
            case 'i':
                rx_centre = atof(options.optarg);
                break;
            case 'j':
                ldpc_en = atoi(options.optarg);
                if ((ldpc_en != 1) && (ldpc_en !=2)) {
                    fprintf(stderr, "--ldpc 1  (224,112) code used for 700D\n");
                    fprintf(stderr, "--ldpc 2  (504,396) code used for 2020\n");
                    exit(1);
                }
                break;
            case 'k':
                txbpf_en = 1;
                break;
            case 'l':
                use_text = 1;
                break;
            case 'p':
                data_bits_per_frame = atoi(options.optarg);
                break;
            case 'q':
                dpsk = 1;
                break;
            case 'v':
                verbose = atoi(options.optarg);
                if (verbose < 0 || verbose > 3)
                    verbose = 0;
        }
    }

    /* Print remaining arguments to give user a hint */

    char *arg;

    while ((arg = optparse_arg(&options)))
        fprintf(stderr, "%s\n", arg);

    if (input_specified) {
        if ((fin = fopen(fin_name, "rb")) == NULL) {
            fprintf(stderr, "Error opening input bits file: %s\n", fin_name);
            exit(-1);
        }
    }

    if (output_specified) {
        if ((fout = fopen(fout_name, "wb")) == NULL) {
            fprintf(stderr, "Error opening output modem sample file: %s\n", fout_name);
            exit(-1);
        }
    }

    if ((ofdm_config = (struct OFDM_CONFIG *) calloc(1, sizeof (struct OFDM_CONFIG))) == NULL) {
        fprintf(stderr, "Out of Memory\n");
        exit(-1);
    }

    ofdm_config->fs = 8000.0f; /* Sample Frequency */
    ofdm_config->ofdm_timing_mx_thresh = 0.30f;
    ofdm_config->ftwindowwidth = 11;
    ofdm_config->bps = 2; /* Bits per Symbol */
    ofdm_config->txtbits = 4; /* number of auxiliary data bits */
    ofdm_config->ns = ns; /* Number of Symbol frames */

    ofdm_config->tx_centre = tx_centre;
    ofdm_config->rx_centre = rx_centre;
    ofdm_config->nc = nc;
    ofdm_config->tcp = tcp;
    ofdm_config->ts = ts;

    ofdm_config->rs = (1.0f / ts); /* Modulating Symbol Rate */

    struct OFDM *ofdm = ofdm_create(ofdm_config);
    assert(ofdm != NULL);

    free(ofdm_config);

    /* Get a copy of the actual modem config */
    ofdm_config = ofdm_get_config_param();

    ofdm_bitsperframe = ofdm_get_bits_per_frame();
    ofdm_nuwbits = (ofdm_config->ns - 1) * ofdm_config->bps - ofdm_config->txtbits;
    ofdm_ntxtbits = ofdm_config->txtbits;

    /* Set up default LPDC code.  We could add other codes here if we like */

    struct LDPC ldpc;
    
    int Nbitsperframe;
    int coded_bits_per_frame;

    if (ldpc_en) {
        if (ldpc_en == 1)
            set_up_hra_112_112(&ldpc, ofdm_config);
        else
            set_up_hra_504_396(&ldpc, ofdm_config);

        /* here is where we can change data bits per frame to a number smaller than LDPC code input data bits_per_frame */

        if (data_bits_per_frame)
            set_data_bits_per_frame(&ldpc, data_bits_per_frame, ofdm_config->bps);
        
        data_bits_per_frame = ldpc.data_bits_per_frame;
        coded_bits_per_frame = ldpc.coded_bits_per_frame;
 
        assert(data_bits_per_frame <= ldpc.ldpc_data_bits_per_frame);
        assert(coded_bits_per_frame <= ldpc.ldpc_coded_bits_per_frame);
        
        if (verbose > 1) {
            fprintf(stderr, "ldpc_data_bits_per_frame = %d\n", ldpc.ldpc_data_bits_per_frame);
            fprintf(stderr, "ldpc_coded_bits_per_frame  = %d\n", ldpc.ldpc_coded_bits_per_frame);
            fprintf(stderr, "data_bits_per_frame = %d\n", data_bits_per_frame);
            fprintf(stderr, "coded_bits_per_frame  = %d\n", coded_bits_per_frame);
            fprintf(stderr, "ofdm_bits_per_frame  = %d\n", ofdm_bitsperframe);
            fprintf(stderr, "interleave_frames: %d\n", interleave_frames);
        }

        assert((ofdm_nuwbits + ofdm_ntxtbits + coded_bits_per_frame) <= ofdm_bitsperframe); /* sanity check */
        
        Nbitsperframe = interleave_frames*data_bits_per_frame;
    } else {
        /* vanilla uncoded input bits mode */
        Nbitsperframe = ofdm_bitsperframe;
    }

    int Nsamperframe = ofdm_get_samples_per_frame();

    if (verbose) {
        ofdm_set_verbose(ofdm, verbose);
        fprintf(stderr, "Nsamperframe: %d, interleave_frames: %d, Nbitsperframe: %d \n",
                Nsamperframe, interleave_frames, Nbitsperframe);
    }

    uint8_t tx_bits_char[Nbitsperframe];
    short tx_scaled[Nsamperframe];
    uint8_t txt_bits_char[ofdm_ntxtbits * interleave_frames];

    for (i = 0; i < ofdm_ntxtbits * interleave_frames; i++) {
        txt_bits_char[i] = 0;
    }

    if (testframes) {
        Nrows = Nsec * ofdm_config->rs;
        Nframes = floor((Nrows - 1) / ofdm_config->ns);

        if (verbose)
            fprintf(stderr, "Nframes: %d\n", Nframes);
    }

    if (txbpf_en) {
        ofdm_set_tx_bpf(ofdm, 1);
    }
    if (dpsk) {
        ofdm_set_dpsk(ofdm, 1);
    }

    char text_str[] = "cq cq cq hello world\r"; // Add text bits to match other tests
    char *ptr_text = text_str;

    short tx_varicode_bits[VARICODE_MAX_BITS];
    int nvaricode_bits = 0;
    int varicode_bit_index = 0;

    if (verbose > 1) {
	ofdm_print_info(ofdm);
    }

    /* main loop ----------------------------------------------------------------*/

    int frame = 0;

    while (fread(tx_bits_char, sizeof (uint8_t), Nbitsperframe, fin) == Nbitsperframe) {

        if (ldpc_en) {
            /* fancy interleaved LDPC encoded frames ----------------------------*/

            /* optionally overwrite input data with test frame of
               payload data bits known to demodulator */

            if (testframes) {

                if (use_text) {
                    // Get text bits
                    int nspare = ofdm_ntxtbits*interleave_frames;
                    int k;

                    for (k = 0; k < nspare; k++) {
                        if (nvaricode_bits) {
                            txt_bits_char[k] = tx_varicode_bits[varicode_bit_index++];
                            nvaricode_bits--;
                        }

                        if (nvaricode_bits == 0) {
                            /* get new char and encode */
                            char s[2];
                            s[0] = *ptr_text++;

                            if (*ptr_text == 0)
                                ptr_text = &text_str[0];

                            nvaricode_bits = varicode_encode(tx_varicode_bits, s, VARICODE_MAX_BITS, 1, 1);
                            varicode_bit_index = 0;
                        }
                    }
                }

                uint16_t r[data_bits_per_frame];

                ofdm_rand(r, data_bits_per_frame);

                for (j = 0; j < interleave_frames; j++) {
                    for (i = 0; i < data_bits_per_frame; i++) {
                        tx_bits_char[j * data_bits_per_frame + i] = r[i] > 16384;
                    }
                }
            }

            complex float tx_sams[interleave_frames * Nsamperframe];
            ofdm_ldpc_interleave_tx(ofdm, &ldpc, tx_sams, tx_bits_char, txt_bits_char, interleave_frames, ofdm_config);

            for (j = 0; j < interleave_frames; j++) {
                for (i = 0; i < Nsamperframe; i++) {
                    tx_scaled[i] = OFDM_AMP_SCALE * crealf(tx_sams[j * Nsamperframe + i]);
                }

                fwrite(tx_scaled, sizeof (short), Nsamperframe, fout);
            }
        } else {
            /* just modulate uncoded raw bits ------------------------------------*/

            if (testframes) {
                /* build up a test frame consisting of unique word, txt bits, and psuedo-random
                   uncoded payload bits.  The psuedo-random generator is the same as Octave so
                   it can interoperate with ofdm_tx.m/ofdm_rx.m */

                int Npayloadbits = Nbitsperframe - (ofdm_nuwbits + ofdm_ntxtbits);
                uint16_t r[Npayloadbits];
                uint8_t payload_bits[Npayloadbits];

                ofdm_rand(r, Npayloadbits);

                for (i = 0; i < Npayloadbits; i++) {
                    payload_bits[i] = r[i] > 16384;
                }

                uint8_t txt_bits[ofdm_ntxtbits];

                for (i = 0; i < ofdm_ntxtbits; i++) {
                    txt_bits[i] = 0;
                }

                ofdm_assemble_modem_frame(ofdm, tx_bits_char, payload_bits, txt_bits);
            }

            int tx_bits[Nbitsperframe];

            for (i = 0; i < Nbitsperframe; i++) {
                tx_bits[i] = tx_bits_char[i];
            }

	    if (verbose >=3) {
                fprintf(stderr, "\ntx_bits:\n");
                for (i = 0; i < Nbitsperframe; i++) {
                    fprintf(stderr, "  %3d %8d\n", i, tx_bits[i]);
                }
            }

            COMP tx_sams[Nsamperframe];
            ofdm_mod(ofdm, tx_sams, tx_bits);

	    if (verbose >=3) {
                fprintf(stderr, "\ntx_sams:\n");
                for (i = 0; i < Nsamperframe; i++) {
                    fprintf(stderr, "  %3d % f\n", i, (double)tx_sams[i].real);
                }
            }

            /* scale and save to disk as shorts */

            for (i = 0; i < Nsamperframe; i++)
                tx_scaled[i] = tx_sams[i].real * OFDM_AMP_SCALE;

            fwrite(tx_scaled, sizeof (short), Nsamperframe, fout);
        }

        frame++;

        if (testframes && (frame >= Nframes))
            break;
    }

    if (input_specified)
        fclose(fin);

    if (output_specified)
        fclose(fout);

    if (verbose)
        fprintf(stderr, "%d frames processed\n", frame);

    ofdm_destroy(ofdm);

    return 0;
}

