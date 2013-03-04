/*---------------------------------------------------------------------------*\

  FILE........: fec_dec.c
  AUTHOR......: David Rowe
  DATE CREATED: 4 march 2013

  FEC decoder for data from modem containing compressed Codec 2 data
  and FEC.

\*---------------------------------------------------------------------------*/

/*
  Copyright (C) 2013 David Rowe

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

#include "codec2.h"
#include "codec2_fdmdv.h"

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>

int main(int argc, char *argv[])
{
    void          *codec2, *fdmdv;
    FILE          *fin;
    FILE          *fout;
    int            bits_per_input_frame, bytes_per_input_frame;
    unsigned char *packed_input_bits;
    int           *unpacked_input_bits;
    int            bits_per_output_frame, bytes_per_output_frame;
    unsigned char *packed_output_bits;
    int           *unpacked_output_bits;
    int            mode, Nc, bit, byte;
    int            i;

    if (argc < 3) {
	printf("%s InputFromModemWithFECFile OutputToCodec2File\n", argv[0]);
	exit(1);
    }

    if (strcmp(argv[1], "-")  == 0) fin = stdin;
    else if ( (fin = fopen(argv[1],"rb")) == NULL ) {
	fprintf(stderr, "Error opening input file from Demod: %s: %s.\n",
         argv[2], strerror(errno));
	exit(1);
    }

    if (strcmp(argv[2], "-") == 0) fout = stdout;
    else if ( (fout = fopen(argv[2],"wb")) == NULL ) {
	fprintf(stderr, "Error opening output file to Codec : %s: %s.\n",
         argv[3], strerror(errno));
	exit(1);
    }

    /* input parameters and buffers. Note data is split into two 20ms
       frames for transmission over modem. */

    Nc = 20;
    fdmdv = fdmdv_create(Nc);

    bits_per_input_frame = 2*fdmdv_bits_per_frame(fdmdv);
    bytes_per_input_frame = bits_per_input_frame / 8;
    assert((bits_per_input_frame % 8) == 0); /* make sure integer number of bytes per frame */

    packed_input_bits = (unsigned char*)malloc(bytes_per_input_frame*sizeof(char));
    assert(packed_input_bits != NULL);
    unpacked_input_bits = (int*)malloc(bits_per_input_frame*sizeof(int));
    assert(unpacked_input_bits != NULL);

    /* 
       Output parameters and buffers.
    */

    mode = CODEC2_MODE_1400;
    codec2 = codec2_create(mode);

    bits_per_output_frame = codec2_bits_per_frame(codec2);
    bytes_per_output_frame = bits_per_output_frame/8;
    assert((bits_per_output_frame % 8) == 0); /* make sure integer number of bytes per frame */

    packed_output_bits = (unsigned char*)malloc(bytes_per_output_frame*sizeof(char));
    assert(packed_output_bits != NULL);
    unpacked_output_bits = (int*)malloc(bits_per_output_frame*sizeof(int));
    assert(unpacked_output_bits != NULL);
    
    fprintf(stderr, "input bits: %d  input_bytes: %d  output_bits: %d  output_bytes: %d\n",
            bits_per_input_frame,  bytes_per_input_frame, bits_per_output_frame,  bytes_per_output_frame);

    /* main loop */

    while(fread(packed_input_bits, sizeof(char), bytes_per_input_frame, fin) == (size_t)bytes_per_input_frame) {

	/* unpack bits, MSB first */

	bit = 7; byte = 0;
	for(i=0; i<bits_per_input_frame; i++) {
	    unpacked_input_bits[i] = (packed_input_bits[byte] >> bit) & 0x1;
	    bit--;
	    if (bit < 0) {
		bit = 7;
		byte++;
	    }
	}
	assert(byte == bytes_per_input_frame);

        for(i=0; i<bits_per_output_frame; i++)
            unpacked_output_bits[i] = unpacked_input_bits[i];

        /* pack bits, MSB first  */

        bit = 7; byte = 0;
        memset(packed_output_bits, 0, bytes_per_output_frame);
        for(i=0; i<bits_per_output_frame; i++) {
            packed_output_bits[byte] |= (unpacked_output_bits[i] << bit);
            bit--;
            if (bit < 0) {
                bit = 7;
                byte++;
            }
        }
        assert(byte == bytes_per_output_frame);
        
 	fwrite(packed_output_bits, sizeof(char), bytes_per_output_frame, fout);

	/* if this is in a pipeline, we probably don't want the usual
           buffering to occur */

        if (fout == stdout) fflush(stdout);
        if (fin == stdin) fflush(stdin);         
    }

    codec2_destroy(codec2);
    fdmdv_destroy(fdmdv);

    free(packed_input_bits);
    free(unpacked_input_bits);
    free(packed_output_bits);
    free(unpacked_output_bits);

    fclose(fin);
    fclose(fout);

    return 0;
}
