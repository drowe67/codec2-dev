/*---------------------------------------------------------------------------*\
                                                                             
  FILE........: fdmdv_demod.c
  AUTHOR......: David Rowe  
  DATE CREATED: April 30 2012
                                                                             
  Given an input raw file (8kHz, 16 bit shorts) of FDMDV modem samples
  outputs a file of bits.  The output file is assumed to be arranged
  as codec frames of 56 bits (7 bytes) which are received as two 28
  bit modem frames.
                                                                             
\*---------------------------------------------------------------------------*/


/*
  Copyright (C) 2012 David Rowe

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

#include "fdmdv.h"
#include "octave.h"

#define BITS_PER_CODEC_FRAME (2*FDMDV_BITS_PER_FRAME)
#define BYTES_PER_CODEC_FRAME (BITS_PER_CODEC_FRAME/8)

/* lof of information we want to dump to Octave */

#define MAX_FRAMES 50*60 /* 1 minute at 50 symbols/s */

int main(int argc, char *argv[])
{
    FILE         *fin, *fout;
    struct FDMDV *fdmdv;
    char          packed_bits[BYTES_PER_CODEC_FRAME];
    int           rx_bits[FDMDV_BITS_PER_FRAME];
    int           codec_bits[2*FDMDV_BITS_PER_FRAME];
    float         rx_fdm[FDMDV_SAMPLES_PER_FRAME];
    short         rx_fdm_scaled[FDMDV_SAMPLES_PER_FRAME];
    int           i, bit, byte, c;
    int           nin;
    int           sync_bit;
    int           state, next_state;

    int           frames;
    FILE         *foct = NULL;
    struct FDMDV_STATS stats;
    COMP          rx_symbols_log[FDMDV_NSYM][MAX_FRAMES];
    int           coarse_fine_log[MAX_FRAMES];
    float         rx_timing_log[MAX_FRAMES];
    float         foff_log[MAX_FRAMES];

    if (argc < 3) {
	printf("usage: %s InputModemRawFile OutputBitFile [OctaveDumpFile]\n", argv[0]);
	printf("e.g    %s hts1a_fdmdv.raw hts1a.c2\n", argv[0]);
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
	fprintf(stderr, "Error opening output bit file: %s: %s.\n",
         argv[2], strerror(errno));
	exit(1);
    }

    fdmdv = fdmdv_create();
    frames = 0;
    state = 0;
    nin = FDMDV_SAMPLES_PER_FRAME;

    while(fread(rx_fdm_scaled, sizeof(short), nin, fin) == nin)
    {
	for(i=0; i<FDMDV_SAMPLES_PER_FRAME; i++)
	    rx_fdm[i] = rx_fdm_scaled[i]/FDMDV_SCALE;
	fdmdv_demod(fdmdv, rx_bits, &sync_bit, rx_fdm, &nin);

	/* log data for optional Octave dump */

	if (frames < MAX_FRAMES) {
	    fdmdv_get_demod_stats(fdmdv, &stats);
	    for(c=0; c<FDMDV_NSYM; c++)
		rx_symbols_log[c][frames] = stats.rx_symbols[c];
	    foff_log[frames] = stats.foff;
	    rx_timing_log[frames] = stats.rx_timing;
	    coarse_fine_log[frames] = stats.fest_coarse_fine;
	    frames++;
	}
	else
	    printf("MAX_FRAMES exceed in Octave log, log truncated\n");

	/* state machine to output codec bits only if we have a 0,1
	   sync bit sequence */

	next_state = state;
	switch (state) {
	case 0:
	    if (sync_bit == 0) {
		next_state = 1;
		memcpy(codec_bits, rx_bits, FDMDV_BITS_PER_FRAME*sizeof(int));
	    }
	    else
		next_state = 0;
	    break;
	case 1:
	    if (sync_bit == 1) {
		memcpy(&codec_bits[FDMDV_BITS_PER_FRAME], rx_bits, FDMDV_BITS_PER_FRAME*sizeof(int));

		/* pack bits, MSB received first  */

		bit = 7; byte = 0;
		memset(packed_bits, 0, BYTES_PER_CODEC_FRAME);
		for(i=0; i<BITS_PER_CODEC_FRAME; i++) {
		    packed_bits[byte] |= (codec_bits[i] << bit);
		    bit--;
		    if (bit < 0) {
			bit = 7;
			byte++;
		    }
		}
		assert(byte == BYTES_PER_CODEC_FRAME);

		fwrite(packed_bits, sizeof(char), BYTES_PER_CODEC_FRAME, fout);
	    }
	    next_state = 0;
	    break;
	}	
	state = next_state;

	/* if this is in a pipeline, we probably don't want the usual
	   buffering to occur */

        if (fout == stdout) fflush(stdout);
        if (fin == stdin) fflush(stdin);         
    }

    /* Optional dump to Octave log file */

    if ( strcmp(argv[3],"|") && (foct = fopen(argv[3],"wt")) == NULL ) {
	fprintf(stderr, "Error opening Octave dump file: %s: %s.\n",
		argv[3], strerror(errno));
	exit(1);
    }
    else {
	octave_save_complex(foct, "rx_symbols_log_c", (COMP*)rx_symbols_log, FDMDV_NSYM, MAX_FRAMES, MAX_FRAMES);  
	octave_save_float(foct, "foff_log_c", foff_log, 1, MAX_FRAMES);  
	octave_save_float(foct, "rx_timing_log_c", rx_timing_log, 1, MAX_FRAMES);  
	octave_save_int(foct, "coarse_fine_log_c", coarse_fine_log, 1, MAX_FRAMES);  
	fclose(foct);
    }

    fclose(fin);
    fclose(fout);
    fdmdv_destroy(fdmdv);

    return 0;
}

