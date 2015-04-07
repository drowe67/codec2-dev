/*---------------------------------------------------------------------------*\
                                                                             
  FILE........: cohpsk_demod.c
  AUTHOR......: David Rowe  
  DATE CREATED: April 6 2015
                                                                             
  Given an input file of raw file (8kHz, 16 bit shorts) of COHPSK modem samples,
  outputs a file of bits (note one bit per int, not compressed).
                                                                             
\*---------------------------------------------------------------------------*/

/*
  Copyright (C) 2015 David Rowe

  All rights reserved.

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU Lesser General Public License version 2, as
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

#include "codec2_cohpsk.h"

int main(int argc, char *argv[])
{
    FILE          *fin, *fout;
    struct COHPSK *cohpsk;
    int           rx_bits[COHPSK_BITS_PER_FRAME];
    COMP          rx_fdm[COHPSK_SAMPLES_PER_FRAME];
    short         rx_fdm_scaled[COHPSK_SAMPLES_PER_FRAME];
    int           frames, reliable_sync_bit;
    int           i;

    if (argc < 3) {
	printf("usage: %s InputModemRawFile OutputOneBitPerIntFile\n", argv[0]);
	exit(1);
    }

    if (strcmp(argv[1], "-")  == 0) fin = stdin;
    else if ( (fin = fopen(argv[1],"rb")) == NULL ) {
	fprintf(stderr, "Error opening input modem sample file: %s: %s.\n",
         argv[1], strerror(errno));
	exit(1);
    }

    if (strcmp(argv[2], "-") == 0) fout = stdout;
    else if ( (fout = fopen(argv[2],"wb")) == NULL ) {
	fprintf(stderr, "Error opening output file: %s: %s.\n",
         argv[2], strerror(errno));
	exit(1);
    }

    cohpsk = cohpsk_create();

    frames = 0;

    while(fread(rx_fdm_scaled, sizeof(short), COHPSK_SAMPLES_PER_FRAME, fin) == COHPSK_SAMPLES_PER_FRAME) {
	frames++;

	/* scale and demod */

	for(i=0; i<COHPSK_SAMPLES_PER_FRAME; i++) {
	    rx_fdm[i].real = rx_fdm_scaled[i]/FDMDV_SCALE;
            rx_fdm[i].imag = 0.0;
        }

	cohpsk_demod(cohpsk, rx_bits, &reliable_sync_bit, rx_fdm);

 	fwrite(rx_bits, sizeof(int), COHPSK_BITS_PER_FRAME, fout);

	/* if this is in a pipeline, we probably don't want the usual
	   buffering to occur */

        if (fout == stdout) fflush(stdout);
        if (fin == stdin) fflush(stdin); 
    }

    fclose(fin);
    fclose(fout);
    cohpsk_destroy(cohpsk);

    return 0;
}
