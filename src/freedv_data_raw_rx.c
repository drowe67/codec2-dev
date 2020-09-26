/*---------------------------------------------------------------------------*\

  FILE........: freedv_data_raw_rx.c
  AUTHOR......: David Rowe
  DATE CREATED: May 2020

  Demonstrates receiving frames of raw data bytes (instead of
  compressed speech) using the FreeDV API.

\*---------------------------------------------------------------------------*/

/*
  Copyright (C) 2020 David Rowe

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
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>

#include "freedv_api.h"

int main(int argc, char *argv[]) {
    FILE                      *fin, *fout;
    struct freedv             *freedv;
    int                        nin, nbytes, nbytes_total = 0, frame = 0;
    int                        mode;
    int                        verbose, use_testframes;
    int                        i;
    
    if (argc < 4) {
        char f2020[80] = {0};
        #ifdef __LPCNET__
        sprintf(f2020,"|2020");
        #endif     
	printf("usage: %s 700C|700D|800XA|FSK_LDPC%s InputModemSpeechFile BinaryDataFile\n"
               " [-v] \n", argv[0],f2020);
	printf("e.g    %s 700D dataBytes_700d.raw dataBytes_rx.bin\n", argv[0]);
	exit(1);
    }

    mode = -1;
    if (!strcmp(argv[1],"700C")) mode = FREEDV_MODE_700C;
    if (!strcmp(argv[1],"700D")) mode = FREEDV_MODE_700D;
    #ifdef __LPCNET__
    if (!strcmp(argv[1],"2020")) mode = FREEDV_MODE_2020;
    #endif
    if (!strcmp(argv[1],"FSK_LDPC")) mode = FREEDV_MODE_FSK_LDPC;
    if (mode == -1) {
        fprintf(stderr, "Error in mode: %s\n", argv[1]);
        exit(1);
    }

    if (strcmp(argv[2], "-")  == 0) fin = stdin;
    else if ( (fin = fopen(argv[2],"rb")) == NULL ) {
	fprintf(stderr, "Error opening input raw modem sample file: %s: %s.\n",
         argv[2], strerror(errno));
	exit(1);
    }

    if (strcmp(argv[3], "-") == 0) fout = stdout;
    else if ( (fout = fopen(argv[3],"wb")) == NULL ) {
	fprintf(stderr, "Error opening output speech sample file: %s: %s.\n",
         argv[3], strerror(errno));
	exit(1);
    }

    verbose = 0; use_testframes = 0;
    
    if (argc > 4) {
        for (i = 4; i < argc; i++) {
            if (strcmp(argv[i], "-v") == 0) verbose = 1;
            else if (strcmp(argv[i], "--testframes") == 0) use_testframes = 1;
            else if (strcmp(argv[i], "-vv") == 0) verbose = 2;
            else {
                fprintf(stderr, "unkown option: %s\n", argv[i]);
                exit(1);
            }
        }
    }

    freedv = freedv_open(mode);
    assert(freedv != NULL);
    freedv_set_verbose(freedv, verbose);
    freedv_set_test_frames(freedv, use_testframes);

    /* for streaming bytes it's much easier use the modes that have a multiple of 8 payload bits/frame */
    assert((freedv_get_bits_per_modem_frame(freedv) % 8) == 0);
    int bytes_per_modem_frame = freedv_get_bits_per_modem_frame(freedv)/8;
    fprintf(stderr, "bytes_per_modem_frame: %d\n", bytes_per_modem_frame);
    uint8_t bytes_out[bytes_per_modem_frame];
    short  demod_in[freedv_get_n_max_modem_samples(freedv)];

    /* We need to work out how many samples the demod needs on each
       call (nin).  This is used to adjust for differences in the tx
       and rx sample clock frequencies */

    nin = freedv_nin(freedv);
    while(fread(demod_in, sizeof(short), nin, fin) == nin) {
        frame++;
        
        nbytes = freedv_rawdatarx(freedv, bytes_out, demod_in);
        nin = freedv_nin(freedv);

        /* Output data if FEC decoding indicates it has no uncorrected bit errors */
        if (!freedv_get_uncorrected_errors(freedv)) {
            fwrite(bytes_out, sizeof(uint8_t), nbytes, fout);
            nbytes_total += nbytes;
        }
        
        if (verbose == 1) {
            fprintf(stderr, "frame: %d nin: %d sync: %d nbytes: %d uncorrected: %d\n",
                    frame, nin, freedv_get_sync(freedv), nbytes, freedv_get_uncorrected_errors(freedv));
        }

	/* if using pipes we probably don't want the usual buffering */
        if (fout == stdout) fflush(stdout);
        if (fin == stdin) fflush(stdin);
    }

    fclose(fin);
    fclose(fout);
    fprintf(stderr, "frames processed: %d  output bytes: %d\n", frame, nbytes_total);

    /* in testframe mode finish up with some stats */
    
    if (freedv_get_test_frames(freedv)) {
        int Tbits = freedv_get_total_bits(freedv);
        int Terrs = freedv_get_total_bit_errors(freedv);
        float uncoded_ber = (float)Terrs/Tbits;
        fprintf(stderr, "BER......: %5.4f Tbits: %5d Terrs: %5d\n", 
		(double)uncoded_ber, Tbits, Terrs);
        if ((mode == FREEDV_MODE_700D) || (mode == FREEDV_MODE_2020) || (mode == FREEDV_MODE_FSK_LDPC)) {
            int Tbits_coded = freedv_get_total_bits_coded(freedv);
            int Terrs_coded = freedv_get_total_bit_errors_coded(freedv);
            float coded_ber = (float)Terrs_coded/Tbits_coded;
            fprintf(stderr, "Coded BER: %5.4f Tbits: %5d Terrs: %5d\n",
                    (double)coded_ber, Tbits_coded, Terrs_coded);

            /* set return code for Ctest */
            if ((uncoded_ber < 0.1f) && (coded_ber < 0.01f))
                return 0;
            else
                return 1;
        }
    }

    freedv_close(freedv);
    return 0;
}

