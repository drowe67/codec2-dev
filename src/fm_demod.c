/*---------------------------------------------------------------------------*\
                                                                             
  FILE........: fm_demod.c
  AUTHOR......: David Rowe  
  DATE CREATED: Feb 24 2015
                                                                             
  Given an input raw file (44.4 kHz, 16 bit shorts) with a FM signal centered 
  11.1 kHz, outputs a file of demodulated audio samples.
                                                                             
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

#include "codec2_fm.h"
#include "octave.h"

#define N 160

int main(int argc, char *argv[])
{
    FILE         *fin, *fout;
    struct FM    *fm;
    short         buf[N];
    float         rx[N];
    float         rx_out[N];
    int           i;

    if (argc < 2) {
	printf("usage: %s InputFMRawFile OutputSpeechRawFile\n", argv[0]);
	printf("e.g    %s fm.raw fm_demodulated.raw\n", argv[0]);
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
	fprintf(stderr, "Error opening output file: %s: %s.\n",
         argv[2], strerror(errno));
	exit(1);
    }

    fm         = fm_create(N);
    fm->Fs     = 44400.0;
    fm->fm_max = 3000.0;
    fm->fd     = 5000.0;
    fm->fc     = fm->Fs/4;

    while(fread(buf, sizeof(short), N, fin) == N) {
	for(i=0; i<N; i++) {
	    rx[i] = buf[i];
        }
        fm_demod(fm, rx_out, rx);
	for(i=0; i<N; i++) {
	    buf[i] = 16384*rx_out[i];
        }
        fwrite(buf, sizeof(short), N, fout);
    }

    fm_destroy(fm);
    fclose(fin);
    fclose(fout);

    return 0;
}
