/*---------------------------------------------------------------------------*\

  FILE........: freedv_rx_profile.c
  AUTHOR......: David Rowe
  DATE CREATED: 13 August 2014

  Profiling freedv_rx() operation on the STM32F4.

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

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>

#include "stm32f4xx_conf.h"
#include "stm32f4xx.h"
#include "gdb_stdio.h"
#include "freedv_api.h"
#include "machdep.h"

#ifdef __EMBEDDED__
#define printf gdb_stdio_printf
#define fopen gdb_stdio_fopen
#define fclose gdb_stdio_fclose
#define fread gdb_stdio_fread
#define fwrite gdb_stdio_fwrite
#endif

int main(int argc, char *argv[]) {
    struct freedv *f;
    short          inbuf[FREEDV_NSAMPLES], outbuf[FREEDV_NSAMPLES];
    FILE          *fin, *fout;
    int            frame, nin, nout = 0;
    PROFILE_VAR(freedv_start);

    machdep_profile_init();

    f = freedv_open(FREEDV_MODE_1600);

    // Transmit ---------------------------------------------------------------------

    frame = 0;

    fin = fopen("mod.raw", "rb");
    if (fin == NULL) {
        printf("Error opening input file\n");
        exit(1);
    }

    fout = fopen("stm_out.raw", "wb");
    if (fout == NULL) {
        printf("Error opening output file\n");
        exit(1);
    }

    nin = freedv_nin(f);
    while (fread(inbuf, sizeof(short), nin, fin) == nin) {
        PROFILE_SAMPLE(freedv_start);
        nout = freedv_rx(f, outbuf, inbuf);
        nin = freedv_nin(f);
        PROFILE_SAMPLE_AND_LOG2(freedv_start, "  demod");     

        machdep_profile_print_logged_samples();
        fwrite(outbuf, sizeof(short), nout, fout);
        printf("frame: %d\n", ++frame);
    }

    fclose(fin);
    fclose(fout);

    return 0;
}

