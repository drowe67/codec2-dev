/*---------------------------------------------------------------------------*\

  FILE........: ofdm_get_test_bits.c
  AUTHOR......: David Rowe
  DATE CREATED: Mar 2018

  Generate input for the OFDM modem in either coded or uncoded mode.

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
#include "ldpc_codes.h"
#include "interldpc.h"
#include "varicode.h"

#define IS_DIR_SEPARATOR(c)     ((c) == '/')

static const char *progname;

void opt_help() {
    fprintf(stderr, "\nUsage: %s [options]\n\n", progname);
    fprintf(stderr, "  --out     filename  Name of OutputOneCharPerBitFile\n");
    fprintf(stderr, "  --frames  n         Number of frames to output (default 10)\n");
    fprintf(stderr, "  --ldpc              Frame length (112) for LDPC (else 238) for Plain (default Plain)\n");
    fprintf(stderr, "  --verbose           Output variable assigned values to stderr\n\n");

    exit(-1);
}

int main(int argc, char *argv[])
{
    struct OFDM  *ofdm;
    struct LDPC  ldpc;
    FILE         *fout;
    char         *fout_name = NULL;
    int          opt, verbose, Nframes, n;
    int          ldpc_en, frames, output_specified;

    char *pn = argv[0] + strlen (argv[0]);

    while (pn != argv[0] && !IS_DIR_SEPARATOR (pn[-1]))
        --pn;
    
    progname = pn;

    /* Turn off stream buffering */

    setvbuf(stdout, NULL, _IONBF, BUFSIZ);

    fout = stdout;
    output_specified = 0;
    frames = 10;
    ldpc_en = 0;
    verbose = 0;

    struct optparse options;

    struct optparse_long longopts[] = {
        {"out",        'o', OPTPARSE_REQUIRED},
        {"frames",     'n', OPTPARSE_REQUIRED},
        {"ldpc",       'l', OPTPARSE_NONE},
        {"verbose",    'v', OPTPARSE_NONE},
        {0, 0, 0}
    };

    optparse_init(&options, argv);

    while ((opt = optparse_long(&options, longopts, NULL)) != -1) {
        switch (opt) {
            case '?':
                opt_help();
            case 'o':
                fout_name = options.optarg;
                output_specified = 1;
                break;
            case 'n':
                frames = atoi(options.optarg);
                break;
            case 'l':
                ldpc_en = 1;
                break;
            case 'v':
                verbose = 1;
        }
    }

    /* Print remaining arguments to give user a hint */

    char *arg;

    while ((arg = optparse_arg(&options)))
        fprintf(stderr, "%s\n", arg);

    if (output_specified) {
        if ((fout = fopen(fout_name, "wb")) == NULL) {
            fprintf(stderr, "Error opening output bit file: %s\n", fout_name);
            exit(-1);
        }
    }

    Nframes = frames;

    if (verbose)
        fprintf(stderr, "Nframes: %d\n", Nframes);

    ofdm = ofdm_create(NULL);
    assert(ofdm != NULL);

    int ofdm_bitsperpacket = ofdm_get_bits_per_packet(ofdm);
    int Ndatabitsperpacket = ofdm_bitsperpacket - ofdm->nuwbits - ofdm->ntxtbits;
    
    /* Optionally set up default LPDC code */
    if (ldpc_en) {
        fprintf(stderr, "codename: %s\n", ofdm->codename);
        ldpc_codes_setup(&ldpc, ofdm->codename);
        Ndatabitsperpacket = ldpc.ldpc_data_bits_per_frame;
    }

    if (verbose)
        fprintf(stderr, "Ndatabitsperpacket: %d\n", Ndatabitsperpacket);

    fprintf(stderr, "Ndatabitsperpacket = %d\n", Ndatabitsperpacket);
    uint8_t data_bits[Ndatabitsperpacket];
    ofdm_generate_payload_data_bits(data_bits, Ndatabitsperpacket);
    for (n = 0; n<Nframes; n++)
	fwrite(data_bits, sizeof(char), Ndatabitsperpacket, fout);

    if (output_specified)
        fclose(fout);

    ofdm_destroy(ofdm);

    return 0;
}

