/*---------------------------------------------------------------------------*\

  FILE........: tcohpsk.c
  AUTHORS.....: David Rowe & Steve Sampson
  DATE CREATED: June 2017

  Tests for the C version of the OFDM modem.  This program
  outputs a file of Octave vectors that are loaded and automatically
  tested against the Octave version of the modem by the Octave script
  tofdm.m

\*---------------------------------------------------------------------------*/

/*
  Copyright (C) 2017 David Rowe

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
#include <complex.h>

#include "ofdm_internal.h"
#include "codec2_ofdm.h"
#include "octave.h"
#include "test_bits_ofdm.h"

#define FRAMES 1

int main(int argc, char *argv[])
{
    struct OFDM   *ofdm;
    COMP           tx[OFDM_SAMPLESPERFRAME];      /* one frame of tx samples */

    int            tx_bits_log[OFDM_BITSPERFRAME*FRAMES];
    COMP           tx_log[OFDM_SAMPLESPERFRAME*FRAMES];

    FILE          *fout;
    int            f;

    ofdm = ofdm_create(); assert(ofdm != NULL);

    /* Main Loop ---------------------------------------------------------------------*/

    for(f=0; f<FRAMES; f++) {

	/* --------------------------------------------------------*\
	                          Mod
	\*---------------------------------------------------------*/

        /* todo: add a longer sequence of test bits through
           ofdm_get/put test bits functin similat to cohpsk/fdmdv */

        ofdm_mod(ofdm, (COMP*)tx, test_bits_ofdm);

        /* tx vector logging */

	memcpy(&tx_bits_log[OFDM_BITSPERFRAME*f], test_bits_ofdm, sizeof(int)*OFDM_BITSPERFRAME);
	memcpy(&tx_log[OFDM_SAMPLESPERFRAME*f], tx, sizeof(COMP)*OFDM_SAMPLESPERFRAME);
    }

    /* --------------------------------------------------------*\
	                        Demod
    \*---------------------------------------------------------*/

    for(f=0; f<FRAMES; f++) {
        /* todo: run demod and log states as it runs */
    }

    /*---------------------------------------------------------*\
               Dump logs to Octave file for evaluation
                      by tofdm.m Octave script
    \*---------------------------------------------------------*/

    fout = fopen("tofdm_out.txt","wt");
    assert(fout != NULL);
    fprintf(fout, "# Created by tofdm.c\n");
    octave_save_complex(fout, "W_c", (COMP*)ofdm->W, OFDM_NC + 2, OFDM_M, OFDM_M);
    octave_save_int(fout, "tx_bits_log_c", tx_bits_log, 1, OFDM_BITSPERFRAME*FRAMES);
    octave_save_complex(fout, "tx_log_c", (COMP*)tx_log, 1, OFDM_SAMPLESPERFRAME*FRAMES,  OFDM_SAMPLESPERFRAME*FRAMES);
    fclose(fout);

    ofdm_destroy(ofdm);

    return 0;
}

