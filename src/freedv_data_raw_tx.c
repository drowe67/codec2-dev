/*---------------------------------------------------------------------------*\

  FILE........: freedv_data_raw_tx.c
  AUTHOR......: David Rowe
  DATE CREATED: May 2020

  Demonstrates transmitting frames of raw data bytes (instead of
  compressed speech) using the FreeDV API and modems.

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
#include <stdint.h>

#include "freedv_api.h"
#include "ofdm_internal.h"

int main(int argc, char *argv[]) {
    FILE                     *fin, *fout;
    struct freedv            *freedv;
    int                       mode;
    int                       use_clip, use_txbpf, testframes, Nframes = 0;
    int                       i;

    if (argc < 4) {
        char f2020[80] = {0};
        #ifdef __LPCNET__
        sprintf(f2020,"|2020");
        #endif     
        printf("usage: %s 700C|700D|800XA|FSK_LDPC%s InputBinaryDataFile OutputModemRawFile\n"
               "  [--clip 0|1] [--txbpf 0|1] [--testframes Nframes]\n", argv[0], f2020);
        printf("e.g    %s 700D dataBytes.bin dataBytes_700d.raw\n", argv[0]);
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
        fprintf(stderr, "Error opening input fle of bytes: %s: %s.\n", argv[2], strerror(errno));
        exit(1);
    }

    if (strcmp(argv[3], "-") == 0) fout = stdout;
    else if ( (fout = fopen(argv[3],"wb")) == NULL ) {
        fprintf(stderr, "Error opening output modem sample file: %s: %s.\n", argv[3], strerror(errno));
        exit(1);
    }

    use_clip = 0; use_txbpf = 0; testframes = 0;
    
    if (argc > 4) {
        for (i = 4; i < argc; i++) {
            if (strcmp(argv[i], "--clip") == 0) { use_clip = atoi(argv[i+1]); i++; }
            else if (strcmp(argv[i], "--txbpf") == 0) { use_txbpf = atoi(argv[i+1]); i++; }
            else if (strcmp(argv[i], "--testframes") == 0) {
                testframes = 1;
                Nframes = atoi(argv[i+1]); i++;
                fprintf(stderr, "Nframes: %d\n", Nframes);
            }
            else {
                fprintf(stderr, "unkown option: %s\n", argv[i]);
                exit(1);
            }
        }
    }

    freedv = freedv_open(mode);
    assert(freedv != NULL);

    /* these are optional ------------------ */
    freedv_set_clip(freedv, use_clip);
    freedv_set_tx_bpf(freedv, use_txbpf);

    /* for streaming bytes it's much easier to use modes that have a multiple of 8 payload bits/frame */
    int bytes_per_modem_frame = freedv_get_bits_per_modem_frame(freedv)/8;
    fprintf(stderr, "bits_per_modem_frame: %d bytes_per_modem_frame: %d\n", freedv_get_bits_per_modem_frame(freedv), bytes_per_modem_frame);
    assert((freedv_get_bits_per_modem_frame(freedv) % 8) == 0);
    int     n_mod_out = freedv_get_n_nom_modem_samples(freedv);
    uint8_t bytes_in[bytes_per_modem_frame];
    short   mod_out[n_mod_out];

    /* optionally set up a known testframe */
    uint8_t testframe_bytes[bytes_per_modem_frame];
    memset(testframe_bytes, 0, bytes_per_modem_frame);
    if (testframes) {
        int bits_per_frame = freedv_get_bits_per_modem_frame(freedv);
        uint8_t testframe_bits[bits_per_frame];
        ofdm_generate_payload_data_bits(testframe_bits, bits_per_frame);
        /* pack bits, MSB first */
        int bit = 7, byte = 0;
        for(int i=0; i<bits_per_frame; i++) {
            testframe_bytes[byte] |= (testframe_bits[i] << bit);
            bit--;
            if (bit < 0) {
                bit = 7;
                byte++;
            }
        }
    }
    fprintf(stderr, "\n");
    
    /* OK main loop  --------------------------------------- */

    int frames = 0;
    while(fread(bytes_in, sizeof(uint8_t), bytes_per_modem_frame, fin) == bytes_per_modem_frame) {
        if (testframes) {
            memcpy(bytes_in, testframe_bytes, bytes_per_modem_frame);
        }

        freedv_rawdatatx(freedv, mod_out, bytes_in);
        fwrite(mod_out, sizeof(short), n_mod_out, fout);
    
        /* if using pipes we don't want the usual buffering to occur */
        if (fout == stdout) fflush(stdout);
        if (fin == stdin) fflush(stdin);

        frames++;       
        if (testframes && (frames >= Nframes)) {
            goto finished;
        }

    }

 finished:
    /* A few extra output buffers so demod can complete */
    for(int i=0; i< n_mod_out; i++) mod_out[i] = 0;
    fwrite(mod_out, sizeof(short), n_mod_out, fout);
    fwrite(mod_out, sizeof(short), n_mod_out, fout);

    freedv_close(freedv);
    fclose(fin);
    fclose(fout);
    
    return 0;
}

