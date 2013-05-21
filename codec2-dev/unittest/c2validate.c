/*---------------------------------------------------------------------------*\

  FILE........: c2validate.c
  AUTHOR......: David Rowe
  DATE CREATED: 10 April 2013

  Encodes and decodes an array of speech samples using Codec 2 and compares
  it to a previously stored output to validate Codec operation.

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
#ifdef __UNITTEST__
#include "hts1a.h"
#include "hts1a_1300.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef __EMBEDDED__
#include "gdb_stdio.h"
#define fopen gdb_stdio_fopen
#define fwrite gdb_stdio_fwrite
#define fclose gdb_stdio_fclose
#endif

int c2validate(int mode, short input_samples[], short output_samples[], char outfile[], int nsamples)
{
    struct CODEC2 *codec2;
    short         *pinput, *poutput, *outbuf;
    unsigned char *bits;
    int            nsam, nbit;
    int            nframes, i, result, j;
    FILE          *fout;

    codec2 = codec2_create(mode);
    nsam = codec2_samples_per_frame(codec2);
    outbuf = (short*)malloc(nsam*sizeof(short));
    nbit = codec2_bits_per_frame(codec2);
    bits = (unsigned char*)malloc(nbit*sizeof(char));
    if (strlen(outfile))
        fout = fopen(outfile, "wb");
    else
        fout = NULL;

    nframes = nsamples/nsam;
    pinput  = input_samples;
    poutput = output_samples;
    result = 1;

    for(i=0; i<nframes; i++) {
	codec2_encode(codec2, bits, pinput);
	codec2_decode(codec2, outbuf, bits);
        for(j=0; j<nsam; j++) {
            if (outbuf[j] != poutput[j])
                result = 0;
        }
        if (fout != NULL)
            fwrite(outbuf, sizeof(short), nsam, fout);
        pinput += nsam;
        poutput += nsam;
    }

    if (fout != NULL)
        fclose(fout);
    free(outbuf);
    free(bits);
    codec2_destroy(codec2);

    return result;
}

#ifdef __UNITTEST__
int main(int argc, char *argv[])
{
    int ret;

    ret = c2validate(CODEC2_MODE_1300, hts1a_raw, hts1a_1300, "", sizeof(hts1a_raw)/sizeof(short));
    if (ret)
        printf("Pass\n");
    else
        printf("Fail\n");
    return 0;
}
#endif
