/*---------------------------------------------------------------------------*\
                                                                             
  FILE........: tfdmdv.c
  AUTHOR......: David Rowe  
  DATE CREATED: April 16 2012
                                                                             
  Unit tests for FDMDV modem.  Combination of unit tests perfromed
  entirely by this program and comparisons with reference Octave
  version of the modem that require running an Octave script
  ../octave/tfdmdv.m.
                                                                             
\*---------------------------------------------------------------------------*/


/*
  Copyright (C) 2009 David Rowe

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

#include "fdmdv_internal.h"
#include "fdmdv.h"

#define FRAMES 10

void octave_save_int(FILE *f, char name[], int data[], int len);
void octave_save_complex(FILE *f, char name[], COMP data[], int len);

int main(int argc, char *argv[])
{
    struct FDMDV *fdmdv;
    int           tx_bits[FDMDV_BITS_PER_FRAME*FRAMES];
    COMP          tx_symbols[(NC+1)*FRAMES];
    FILE         *fout;
    int           f,i;

    fdmdv = fdmdv_create();

    for(f=0; f<FRAMES; f++) {
	fdmdv_get_test_bits(fdmdv, &tx_bits[FDMDV_BITS_PER_FRAME*f]);
	bits_to_dqpsk_symbols(&tx_symbols[(NC+1)*f], fdmdv->prev_tx_symbols, 
			      &tx_bits[FDMDV_BITS_PER_FRAME*f], &fdmdv->tx_pilot_bit);
	memcpy(fdmdv->prev_tx_symbols, &tx_symbols[(NC+1)*f], (NC+1)*sizeof(COMP));
    }
 
    codec2_destroy(fdmdv);

    /* dump to Octave file for evaluation by Octave script */

    fout = fopen("tfdmdv_out.txt","wt");
    assert(fout != NULL);
    fprintf(fout, "# Created by tfdmdv.c\n");
    octave_save_int(fout, "tx_bits_tfdmdv", tx_bits, FDMDV_BITS_PER_FRAME*FRAMES);
    octave_save_complex(fout, "tx_symbols_tfdmdv", tx_symbols, (NC+1)*FRAMES);  
    fclose(fout);

    return 0;
}

void octave_save_int(FILE *f, char name[], int data[], int len)
{
    int i;

    fprintf(f, "# name: %s\n", name);
    fprintf(f, "# type: matrix\n");
    fprintf(f, "# rows: %d\n", 1);
    fprintf(f, "# columns: %d\n", len);
    
    for(i=0; i<len; i++)
	fprintf(f, " %d", data[i]);

    fprintf(f, "\n\n\n");
}

void octave_save_complex(FILE *f, char name[], COMP data[], int len)
{
    int i;

    fprintf(f, "# name: %s\n", name);
    fprintf(f, "# type: complex matrix\n");
    fprintf(f, "# rows: %d\n", 1);
    fprintf(f, "# columns: %d\n", len);
    
    for(i=0; i<len; i++)
	fprintf(f, " (%f,%f)", data[i].real, data[i].imag);

    fprintf(f, "\n\n\n");
}
