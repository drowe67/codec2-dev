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

#include "codec2_ofdm.h"

#define ASCALE (2E5*1.1491)

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
    FILE          *fin, *fout;
    struct OFDM   *ofdm;
    int           frames;
    int           i;

    if (argc < 3) {
        fprintf(stderr, "\n");
	fprintf(stderr, "usage: %s InputOneCharPerBitFile OutputModemRawFile\n", argv[0]);
        fprintf(stderr, "\n");
	exit(1);
    }

    if (strcmp(argv[1], "-")  == 0) fin = stdin;
    else if ( (fin = fopen(argv[1],"rb")) == NULL ) {
	fprintf(stderr, "Error opening input file: %s: %s.\n",
         argv[1], strerror(errno));
	exit(1);
    }

    if (strcmp(argv[2], "-") == 0) fout = stdout;
    else if ( (fout = fopen(argv[2],"wb")) == NULL ) {
	fprintf(stderr, "Error opening output modem sample file: %s: %s.\n",
         argv[2], strerror(errno));
	exit(1);
    }

    ofdm = ofdm_create(OFDM_CONFIG_700D);
    assert(ofdm != NULL);
    int Nbitsperframe = ofdm_get_bits_per_frame(ofdm);
    int Nsamperframe = ofdm_get_samples_per_frame();

    char  tx_bits_char[Nbitsperframe];
    int   tx_bits[Nbitsperframe];
    COMP  tx_sams[Nsamperframe];
    short tx_scaled[Nsamperframe];
    
    frames = 0;

    while(fread(tx_bits_char, sizeof(char), Nbitsperframe, fin) == Nbitsperframe) {
	frames++;
        
        for(i=0; i<Nbitsperframe; i++)
            tx_bits[i] = tx_bits_char[i];
	ofdm_mod(ofdm, tx_sams, tx_bits);

	/* scale and save to disk as shorts */

	for(i=0; i<Nsamperframe; i++)
	    tx_scaled[i] = ASCALE * tx_sams[i].real;

 	fwrite(tx_scaled, sizeof(short), Nsamperframe, fout);

	/* if this is in a pipeline, we probably don't want the usual
	   buffering to occur */

        if (fout == stdout) fflush(stdout);
        if (fin == stdin) fflush(stdin);
    }

    fclose(fin);
    fclose(fout);
    ofdm_destroy(ofdm);

    return 0;
}
