/*---------------------------------------------------------------------------*\
                                                                             
  FILE........: freedv_tx.c
  AUTHOR......: David Rowe
  DATE CREATED: August 2014
                                                                             
  Demo transmit program for FreeDV API functions.
                                                                     
\*---------------------------------------------------------------------------*/

/*
  Copyright (C) 2014 David Rowe

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
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <errno.h>

#include "freedv_api.h"

int main(int argc, char *argv[]) {
    FILE          *fin, *fout;
    short          speech_in[FREEDV_NSAMPLES];
    short          mod_out[FREEDV_NSAMPLES];
    struct freedv *freedv;

    if (argc < 3) {
	printf("usage: %s InputRawSpeechFile OutputModemRawFile\n", argv[0]);
	printf("e.g    %s hts1a.raw hts1a_fdmdv.raw\n", argv[0]);
	exit(1);
    }

    if (strcmp(argv[1], "-")  == 0) fin = stdin;
    else if ( (fin = fopen(argv[1],"rb")) == NULL ) {
	fprintf(stderr, "Error opening input raw speech sample file: %s: %s.\n",
         argv[1], strerror(errno));
	exit(1);
    }

    if (strcmp(argv[2], "-") == 0) fout = stdout;
    else if ( (fout = fopen(argv[2],"wb")) == NULL ) {
	fprintf(stderr, "Error opening output modem sample file: %s: %s.\n",
         argv[2], strerror(errno));
	exit(1);
    }
    
    freedv = freedv_open(FREEDV_MODE_1600);
    assert(freedv != NULL);

    while(fread(speech_in, sizeof(short), FREEDV_NSAMPLES, fin) == FREEDV_NSAMPLES) {
        freedv_tx(freedv, mod_out, speech_in);
        fwrite(mod_out, sizeof(short), FREEDV_NSAMPLES, fout);
    }

    freedv_close(freedv);
    fclose(fin);
    fclose(fout);

    return 0;
}

