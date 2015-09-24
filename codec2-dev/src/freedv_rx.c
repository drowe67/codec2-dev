/*---------------------------------------------------------------------------*\

  FILE........: freedv_rx.c
  AUTHOR......: David Rowe
  DATE CREATED: August 2014

  Demo receive program for FreeDV API functions, some side information
  written to freedv_rx_log.txt

  Example usage (all one line):

     codec2-dev/build_linux/src$ ./freedv_tx 1600 ../../raw/ve9qrp_10s.raw - |
                                 ./freedv_rx 1600 - - | play -t raw -r 8000 -s -2 -

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
#include <stdlib.h>
#include <stdio.h>

#include "freedv_api.h"

struct my_callback_state {
    FILE *ftxt;
};

void my_put_next_rx_char(void *callback_state, char c) {
    struct my_callback_state* pstate = (struct my_callback_state*)callback_state;
    if (pstate->ftxt != NULL) {
        //fprintf(pstate->ftxt, "%c\n", c);
    }
}

int main(int argc, char *argv[]) {
    FILE                      *fin, *fout, *ftxt;
    short                     *speech_out;
    short                     *demod_in;
    struct freedv             *freedv;
    int                        nin, nout, frame = 0;
    struct my_callback_state   my_cb_state;
    int                        mode;
    int                        sync;
    int                        total_bits;
    int                        total_bit_errors;
    float                      snr_est;
    int                        n_speech_samples;
    int                        n_max_modem_samples;

    if (argc < 4) {
	printf("usage: %s 1600|700|700B InputModemSpeechFile OutputSpeechRawFile [--test_frames]\n", argv[0]);
	printf("e.g    %s 1600 hts1a_fdmdv.raw hts1a_out.raw txtLogFile\n", argv[0]);
	exit(1);
    }

    mode = -1;
    if (!strcmp(argv[1],"1600"))
        mode = FREEDV_MODE_1600;
    if (!strcmp(argv[1],"700"))
        mode = FREEDV_MODE_700;
    if (!strcmp(argv[1],"700B"))
        mode = FREEDV_MODE_700B;
    assert(mode != -1);

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

    freedv = freedv_open(mode);
    assert(freedv != NULL);

    if ( (argc > 4) && (strcmp(argv[4], "--testframes") == 0) ) {
      freedv_set_test_frames(freedv, 1);
    }
    freedv_set_snr_squelch_thresh(freedv, -100.0);
    freedv_set_squelch_en(freedv, 1);

    n_speech_samples = freedv_get_n_speech_samples(freedv);
    n_max_modem_samples = freedv_get_n_max_modem_samples(freedv);
    speech_out = (short*)malloc(sizeof(short)*n_speech_samples);
    assert(speech_out != NULL);
    demod_in = (short*)malloc(sizeof(short)*n_max_modem_samples);
    assert(demod_in != NULL);

    ftxt = fopen("freedv_rx_log.txt","wt");
    assert(ftxt != NULL);
    my_cb_state.ftxt = ftxt;
    freedv_set_callback_txt(freedv, &my_put_next_rx_char, NULL, &my_cb_state);

    /* Note we need to work out how many samples demod needs on each
       call (nin).  This is used to adjust for differences in the tx and rx
       sample clock frequencies.  Note also the number of output
       speech samples is time varying (nout). */

    nin = freedv_nin(freedv);
    while(fread(demod_in, sizeof(short), nin, fin) == nin) {
        frame++;

        nout = freedv_rx(freedv, speech_out, demod_in);
        nin = freedv_nin(freedv);

        fwrite(speech_out, sizeof(short), nout, fout);
        freedv_get_modem_stats(freedv, &sync, &snr_est);
        total_bit_errors = freedv_get_total_bit_errors(freedv);

        /* log some side info to the txt file */

        if (ftxt != NULL) {
            fprintf(ftxt, "frame: %d  demod sync: %d  nin:%d demod snr: %3.2f dB  bit errors: %d\n",
                    frame, sync, nin, snr_est, total_bit_errors);
        }

	/* if this is in a pipeline, we probably don't want the usual
           buffering to occur */

        if (fout == stdout) fflush(stdout);
        if (fin == stdin) fflush(stdin);
    }

    if (freedv_get_test_frames(freedv)) {
        total_bits = freedv_get_total_bits(freedv);
        total_bit_errors = freedv_get_total_bit_errors(freedv);
        fprintf(stderr, "bits: %d errors: %d BER: %3.2f\n", total_bits, total_bit_errors, (float)total_bit_errors/total_bits);
    }

    free(speech_out);
    free(demod_in);
    freedv_close(freedv);
    fclose(fin);
    fclose(fout);

    return 0;
}

