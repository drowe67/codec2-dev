/*---------------------------------------------------------------------------*\

  FILE........: ofdm_gen_test_bits.c
  AUTHOR......: David Rowe
  DATE CREATED: Mar 2018

  A modified version of ofdm_get_test_bits which computes
  the bits using the new ofdm_rand function.
  This version supports either LDPC or Plain format, and 
  inserts both ?? and txt bits, copied from ofdm_mod (which was copied from freedv_tx).

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

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <errno.h>

#include "codec2_ofdm.h"
#include "ofdm_internal.h"
#include "interldpc.h"
#include "varicode.h"

static struct OFDM_CONFIG *ofdm_config;

static int ofdm_bitsperframe;
static int ofdm_nuwbits;
static int ofdm_ntxtbits;

int opt_exists(char *argv[], int argc, char opt[]) {
    int i;
    for (i=0; i<argc; i++) {
        if (strcmp(argv[i], opt) == 0) {
            return i;
        }
    }
    return 0;
}

int main(int argc, char *argv[])
{
    struct OFDM  *ofdm;
    struct LDPC   ldpc;
    FILE         *fout;
    int           Nframes, i, n;
    int           ldpc_en;

    if (argc < 3) {
	fprintf(stderr, "usage: %s OutputOneCharPerBitFile numFrames [--ldpc]\n", argv[0]);
	fprintf(stderr, " --ldpc         Length (238) for LDPC (else plain, 224)\n");
	exit(1);
    }

    if (strcmp(argv[1], "-") == 0)
        fout = stdout;
    else if ( (fout = fopen(argv[1],"wb")) == NULL ) {
	fprintf(stderr, "Error opening output file: %s: %s.\n",
         argv[1], strerror(errno));
	exit(1);
    }

    Nframes = atoi(argv[2]);
    fprintf(stderr, "Nframes: %d\n", Nframes);

    ldpc_en = 0;
    if (opt_exists(argv, argc, "--ldpc")) ldpc_en = 1;


    // Build OFDM and LDPC to get sizes
    if ((ofdm_config = (struct OFDM_CONFIG *) calloc(1, sizeof (struct OFDM_CONFIG))) == NULL) {
	printf("Out of Memory");
	exit(1);
    }

    ofdm = ofdm_create(ofdm_config);
    assert(ofdm != NULL);

    free(ofdm_config);

    /* Get a copy of the actual modem config */
    ofdm_config = ofdm_get_config_param();

    ofdm_bitsperframe = ofdm_get_bits_per_frame();
    ofdm_nuwbits = (ofdm_config->ns - 1) * ofdm_config->bps - ofdm_config->txtbits;
    ofdm_ntxtbits = ofdm_config->txtbits;

    /* Set up default LPDC code.  We could add other codes here if we like */
    set_up_hra_112_112(&ldpc, ofdm_config);

    int data_bits_per_frame = ldpc.data_bits_per_frame;

    int Nbitsperframe;
    if (ldpc_en) { Nbitsperframe = data_bits_per_frame; }
    else { Nbitsperframe = ofdm_bitsperframe; }
    fprintf(stderr, "Nbitsperframe: %d\n", Nbitsperframe);

    uint8_t tx_bits[Nbitsperframe];
    uint8_t txt_bits[ofdm_ntxtbits];

    for(i=0; i< ofdm_ntxtbits; i++) {
        txt_bits[i] = 0;
    }


    // Add text bits to match other tests
    char   text_str[] = "cq cq cq hello world\r";
    char  *ptr_text = &text_str[0];
    short  tx_varicode_bits[VARICODE_MAX_BITS];
    int    nvaricode_bits = 0;
    int    varicode_bit_index = 0;


    for(n=0; n<Nframes; n++) {

        if (ldpc_en) { /* fancy interleaved LDPC encoded frames */

            /* build up a test frame consisting of unique word, txt bits, and 
	    do-random uncoded payload bits.  The psuedo-random generator 
	    he same as Octave so it can interoperate with ofdm_tx.m/ofdm_rx.m */

	    // Get text bits
            int nspare = ofdm_ntxtbits;
	    int k;

            for(k=0; k<nspare; k++) {
                if (nvaricode_bits) {
                    txt_bits[k] = tx_varicode_bits[varicode_bit_index++];
                    nvaricode_bits--;
                }
                if (nvaricode_bits == 0) {
                    /* get new char and encode */
                    char s[2];
                    s[0] = *ptr_text++;
                    if (*ptr_text == 0) ptr_text = &text_str[0];
                    nvaricode_bits = varicode_encode(tx_varicode_bits, s, VARICODE_MAX_BITS, 1, 1);
                    varicode_bit_index = 0;
                }
	        fprintf(stderr, "txt_bits[%d] = %d\n", k, txt_bits[k]);
            }

            uint16_t r[data_bits_per_frame];

            ofdm_rand(r, data_bits_per_frame);

            for(i=0; i<data_bits_per_frame; i++) {
                tx_bits[i] = r[i]>16384;
            }


        } else { // (!ldpc_en)

            int Npayloadbits = Nbitsperframe-(ofdm_nuwbits+ofdm_ntxtbits);
            uint16_t r[Npayloadbits];
            uint8_t  payload_bits[Npayloadbits];

            ofdm_rand(r, Npayloadbits);

            for(i=0; i<Npayloadbits; i++) {
                payload_bits[i] = r[i] > 16384;
            }

            ofdm_assemble_modem_frame(tx_bits, payload_bits, txt_bits);
	}

	fwrite(tx_bits, sizeof(char), Nbitsperframe, fout);

	/* if this is in a pipeline, we probably don't want the usual
	   buffering to occur */

        if (fout == stdout)
            fflush(stdout);
    }

    fclose(fout);

    ofdm_destroy(ofdm);

    return 0;
}

