/*---------------------------------------------------------------------------*\

  FILE........: freedv_datac0c1_rx.c
  AUTHOR......: David Rowe
  DATE CREATED: Dec 2021

  Demonstrates receiving frames of raw data bytes using the FreeDV API.  Two
  parallel receivers are running, so we can receive either DATAC0 or DATAC1
  frames.

  See freedv_datac0c1_tx.c for instructions.
  
\*---------------------------------------------------------------------------*/

/*
  Copyright (C) 2021 David Rowe

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
#include <stdint.h>
#include <string.h>

#include "freedv_api.h"

#define NBUF 160

int main(int argc, char *argv[]) {

    // set up DATAC0 Rx
    struct freedv *freedv_c0 = freedv_open(FREEDV_MODE_DATAC0);
    assert(freedv_c0 != NULL);
    freedv_set_frames_per_burst(freedv_c0, 1);
    freedv_set_verbose(freedv_c0, 0);
    int bytes_per_modem_frame_c0 = freedv_get_bits_per_modem_frame(freedv_c0)/8;
    uint8_t bytes_out_c0[bytes_per_modem_frame_c0];
    short  demod_in_c0[freedv_get_n_max_modem_samples(freedv_c0)];
    
    // set up DATAC1 Rx
    struct freedv *freedv_c1 = freedv_open(FREEDV_MODE_DATAC1);
    assert(freedv_c1 != NULL);
    freedv_set_frames_per_burst(freedv_c1, 1);
    freedv_set_verbose(freedv_c1, 0);
    int bytes_per_modem_frame_c1 = freedv_get_bits_per_modem_frame(freedv_c1)/8;
    uint8_t bytes_out_c1[bytes_per_modem_frame_c1];
    short  demod_in_c1[freedv_get_n_max_modem_samples(freedv_c1)];

    size_t n_c0 = 0;
    size_t n_c1 = 0;
    size_t c0_frames = 0;
    size_t c1_frames = 0;
    short buf[NBUF];
    size_t nin;
    
    while(fread(buf, sizeof(short), NBUF, stdin) == NBUF) {

        // NBUF new samples into DATAC1 Rx
        memcpy(&demod_in_c1[n_c1], buf, sizeof(short)*NBUF);
        n_c1 += NBUF; assert(n_c1 <= freedv_get_n_max_modem_samples(freedv_c1));
        nin = freedv_nin(freedv_c1);
        if (n_c1 > nin) {
            size_t nbytes_out = freedv_rawdatarx(freedv_c1, bytes_out_c1, demod_in_c1);
            if (nbytes_out) {
                fprintf(stderr, "DATAC1 frame received!\n");
                c1_frames++;
            }
            // nin samples were read
            n_c1 -= nin; assert(n_c1 >= 0);
            memmove(demod_in_c1, &demod_in_c1[nin], sizeof(short)*n_c1);
        }
    }

    fprintf(stderr, "DATAC0 Frames: %ld\n", c0_frames);
    fprintf(stderr, "DATAC1 Frames: %ld\n", c1_frames);
    freedv_close(freedv_c0);
    freedv_close(freedv_c1);

    return 0;
}
