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

#define NFRAMES 3

int main(int argc, char *argv[])
{
    struct OFDM   *ofdm;
    COMP           tx[OFDM_SAMPLESPERFRAME];      /* one frame of tx samples */
    int            rx_bits[OFDM_BITSPERFRAME];    /* one frame of rx bits    */

    /* log arrays */

    int            tx_bits_log[OFDM_BITSPERFRAME*NFRAMES];
    COMP           tx_log[OFDM_SAMPLESPERFRAME*NFRAMES];
    COMP           rxbuf_in_log[OFDM_SAMPLESPERFRAME*NFRAMES];
    COMP           rxbuf_log[OFDM_RXBUF*NFRAMES];

    FILE          *fout;
    int            f,i;

    ofdm = ofdm_create(); assert(ofdm != NULL);

    /* Main Loop ---------------------------------------------------------------------*/

    for(f=0; f<NFRAMES; f++) {

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

    COMP *rx = tx_log;

    /* Init rx with ideal timing so we can test with timing estimation disabled */

    int  Nsam = OFDM_SAMPLESPERFRAME*NFRAMES;
    int  prx = 0;
    int  nin = OFDM_SAMPLESPERFRAME + 2*(OFDM_M+OFDM_NCP);

    int  lnew;
    COMP rxbuf_in[OFDM_SAMPLESPERFRAME];

    for (i=0; i<nin; i++,prx++) {
         ofdm->rxbuf[OFDM_RXBUF-nin+i] = rx[prx].real + I*rx[prx].imag;
    }

    int nin_tot = 0;

    for(f=0; f<NFRAMES; f++) {
        /* For initial testng, timing est is off, so nin is always
           fixed.  TODO: we need a constant for rxbuf_in[] size that
           is the maximum possible nin */

        nin = ofdm->nin;
        assert(nin == OFDM_SAMPLESPERFRAME);

        /* Insert samples at end of buffer, set to zero if no samples
           available to disable phase estimation on future pilots on
           last frame of simulation. */

        if ((Nsam-prx) < nin) {
            lnew = Nsam-prx;
        } else {
            lnew = nin;
        }
        //printf("nin: %d prx: %d lnew: %d\n", nin, prx, lnew);
        for(i=0; i<nin; i++) {
            rxbuf_in[i].real = 0.0;
            rxbuf_in[i].imag = 0.0;
        }

        if (lnew) {
            for(i=0; i<lnew; i++, prx++) {
                rxbuf_in[i] = rx[prx];
            }
        }
        assert(prx <= OFDM_SAMPLESPERFRAME*NFRAMES);
#ifdef T

        //ofdm_demod(ofdm, rx_bits, rxbuf_in);

        /* rx vector logging -----------------------------------*/

#endif
        assert(nin_tot < OFDM_SAMPLESPERFRAME*NFRAMES);
	memcpy(&rxbuf_in_log[nin_tot], rxbuf_in, sizeof(COMP)*nin);
        nin_tot += nin;

        for(i=0; i<OFDM_RXBUF; i++) {
            rxbuf_log[OFDM_RXBUF*f+i].real = crealf(ofdm->rxbuf[i]);
            rxbuf_log[OFDM_RXBUF*f+i].imag = cimagf(ofdm->rxbuf[i]);
       }
    }

    /*---------------------------------------------------------*\
               Dump logs to Octave file for evaluation
                      by tofdm.m Octave script
    \*---------------------------------------------------------*/

    fout = fopen("tofdm_out.txt","wt");
    assert(fout != NULL);
    fprintf(fout, "# Created by tofdm.c\n");
    octave_save_complex(fout, "W_c", (COMP*)ofdm->W, OFDM_NC + 2, OFDM_M, OFDM_M);
    octave_save_int(fout, "tx_bits_log_c", tx_bits_log, 1, OFDM_BITSPERFRAME*NFRAMES);
    octave_save_complex(fout, "tx_log_c", (COMP*)tx_log, 1, OFDM_SAMPLESPERFRAME*NFRAMES,  OFDM_SAMPLESPERFRAME*NFRAMES);
    octave_save_complex(fout, "rxbuf_in_log_c", (COMP*)rxbuf_in_log, 1, nin_tot, nin_tot);
    octave_save_complex(fout, "rxbuf_log_c", (COMP*)rxbuf_log, 1, OFDM_RXBUF*NFRAMES,  OFDM_RXBUF*NFRAMES);
    fclose(fout);

    ofdm_destroy(ofdm);

    return 0;
}

