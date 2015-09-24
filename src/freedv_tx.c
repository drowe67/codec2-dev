/*---------------------------------------------------------------------------*\

  FILE........: freedv_tx.c
  AUTHOR......: David Rowe
  DATE CREATED: August 2014

  Demo transmit program for FreeDV API functions.

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

#include "freedv_api.h"

struct my_callback_state {
    char  tx_str[80];
    char *ptx_str;
};

char my_get_next_tx_char(void *callback_state) {
    struct my_callback_state* pstate = (struct my_callback_state*)callback_state;
    char  c = *pstate->ptx_str++;

    if (*pstate->ptx_str == 0) {
        pstate->ptx_str = pstate->tx_str;
    }

    return c;
}

int main(int argc, char *argv[]) {
    FILE                     *fin, *fout;
    short                    *speech_in;
    short                    *mod_out;
    struct freedv            *freedv;
    struct my_callback_state  my_cb_state;
    int                       mode;
    int                       n_speech_samples;
    int                       n_nom_modem_samples;

    if (argc < 4) {
	printf("usage: %s 1600|700|700B InputRawSpeechFile OutputModemRawFile [--testframes]\n", argv[0]);
	printf("e.g    %s 1600 hts1a.raw hts1a_fdmdv.raw\n", argv[0]);
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
	fprintf(stderr, "Error opening input raw speech sample file: %s: %s.\n",
         argv[2], strerror(errno));
	exit(1);
    }

    if (strcmp(argv[3], "-") == 0) fout = stdout;
    else if ( (fout = fopen(argv[3],"wb")) == NULL ) {
	fprintf(stderr, "Error opening output modem sample file: %s: %s.\n",
         argv[3], strerror(errno));
	exit(1);
    }

    freedv = freedv_open(mode);
    assert(freedv != NULL);

    if ((argc > 4) && (strcmp(argv[4], "--testframes") == 0)) {
		freedv_set_test_frames(freedv, 1);
    }
    freedv_set_snr_squelch_thresh(freedv, -100.0);
    freedv_set_squelch_en(freedv, 1);

    n_speech_samples = freedv_get_n_speech_samples(freedv);
    n_nom_modem_samples = freedv_get_n_nom_modem_samples(freedv);
    speech_in = (short*)malloc(sizeof(short)*n_speech_samples);
    assert(speech_in != NULL);
    mod_out = (short*)malloc(sizeof(short)*n_nom_modem_samples);
    assert(mod_out != NULL);

    /* set up callback for txt msg chars */

    sprintf(my_cb_state.tx_str, "cq cq cq hello world\n");
    my_cb_state.ptx_str = my_cb_state.tx_str;
    freedv_set_callback_txt(freedv, NULL, &my_get_next_tx_char, &my_cb_state);

    /* OK main loop */

    while(fread(speech_in, sizeof(short), n_speech_samples, fin) == n_speech_samples) {
        freedv_tx(freedv, mod_out, speech_in);
        fwrite(mod_out, sizeof(short), n_nom_modem_samples, fout);

	/* if this is in a pipeline, we probably don't want the usual
           buffering to occur */

        if (fout == stdout) fflush(stdout);
        if (fin == stdin) fflush(stdin);
    }

    free(speech_in);
    free(mod_out);
    freedv_close(freedv);
    fclose(fin);
    fclose(fout);

    return 0;
}

