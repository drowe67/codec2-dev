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

int main(int argc, char *argv[])
{
    struct OFDM  *ofdm;
    FILE         *fout;
    int           Nsec, Nrows, Nframes, i,n;

    if (argc < 2) {
	printf("usage: %s OutputOneCharPerBitFile numSecs\n", argv[0]);
	exit(1);
    }

    if (strcmp(argv[1], "-") == 0) fout = stdout;
    else if ( (fout = fopen(argv[1],"wb")) == NULL ) {
	fprintf(stderr, "Error opening output file: %s: %s.\n",
         argv[1], strerror(errno));
	exit(1);
    }

    ofdm = ofdm_create(OFDM_CONFIG_700D);
    assert(ofdm != NULL);
    int Nbitsperframe = ofdm_get_bits_per_frame(ofdm);
    char  tx_bits_char[Nbitsperframe];
    ofdm_destroy(ofdm);

    for(i=0; i<Nbitsperframe; i++) {
        tx_bits_char[i] = test_bits_ofdm[i];
    }
    
    Nsec = atoi(argv[2]);
    Nrows = Nsec*OFDM_RS;
    Nframes = floor((Nrows-1)/OFDM_NS);

    for(n=0; n<Nframes; n++) {

	fwrite(tx_bits_char, sizeof(char), Nbitsperframe, fout);

	/* if this is in a pipeline, we probably don't want the usual
	   buffering to occur */

        if (fout == stdout) fflush(stdout);
    }

    fclose(fout);

    return 0;
}
