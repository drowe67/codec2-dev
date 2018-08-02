/*---------------------------------------------------------------------------*\

  FILE........: ofdm_get_test_bits.c
  AUTHOR......: David Rowe
  DATE CREATED: Mar 2018

  Generates frames of test bits, useful for input to ofdm_mod.

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
#include "test_bits_ofdm.h"

static struct OFDM_CONFIG *ofdm_config;

static int ofdm_bitsperframe;
static int ofdm_nuwbits;
static int ofdm_ntxtbits;

int main(int argc, char *argv[])
{
    struct OFDM  *ofdm;
    FILE         *fout;
    int           Nframes, i, n;

    if (argc < 2) {
	printf("usage: %s OutputOneCharPerBitFile [-f] numSecsorFrames\n", argv[0]);
	exit(1);
    }

    if (strcmp(argv[1], "-") == 0)
        fout = stdout;
    else if ( (fout = fopen(argv[1],"wb")) == NULL ) {
	fprintf(stderr, "Error opening output file: %s: %s.\n",
         argv[1], strerror(errno));
	exit(1);
    }

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

    char  tx_bits_char[ofdm_bitsperframe];


    for(i=0; i<ofdm_bitsperframe; i++) {
        tx_bits_char[i] = test_bits_ofdm[i];
    }

    if (strcmp(argv[2], "-f") == 0) {
        Nframes = atoi(argv[3]);
    } else {
        int Nsec = atoi(argv[2]);
        int Nrows = (int)(Nsec * ofdm_config->rs);
        Nframes = floorf((Nrows-1)/ofdm_config->ns);
        fprintf(stderr, "ofdm_bitsperframe: %d Nsec: %d Nrows: %d Nframes: %d\n", ofdm_bitsperframe, Nsec, Nrows, Nframes);
    }


    for(n=0; n<Nframes; n++) {

	fwrite(tx_bits_char, sizeof(char), ofdm_bitsperframe, fout);

	/* if this is in a pipeline, we probably don't want the usual
	   buffering to occur */

        if (fout == stdout)
            fflush(stdout);
    }

    fclose(fout);

    ofdm_destroy(ofdm);

    return 0;
}

