/*---------------------------------------------------------------------------*\

  FILE........: freedv_data_raw_tx.c
  AUTHOR......: David Rowe
  DATE CREATED: May 2020

  Demonstrates transmitting frames of raw data bytes (instead of
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
#include <stdint.h>
#include <getopt.h>

#include "freedv_api.h"
#include "fsk.h"
#include "ofdm_internal.h"

void comp_to_short(short mod_out_short[], COMP mod_out_comp[], int n_mod_out) {
    for(int i=0; i<n_mod_out; i++) {
        mod_out_short[2*i] = (short)(mod_out_comp[i].real);
        mod_out_short[2*i+1] = (short)(mod_out_comp[i].imag);
    }            
}

int main(int argc, char *argv[]) {
    FILE                     *fin, *fout;
    struct freedv_advanced    adv = {0,2,100,8000,1000,200, "H_256_512_4"};
    struct freedv            *freedv;
    int                       mode;
    int                       use_clip, use_txbpf, testframes, Nframes = 0;
    int                       use_complex = 0;
    float                     amp = FSK_SCALE;
    int                       shorts_per_sample = 1;
    int                       Nbursts = 1, sequence_numbers = 0;
    
    char f2020[80] = {0};
    if (argc < 4) {
    helpmsg:
        #ifdef __LPCNET__
        sprintf(f2020,"|2020");
        #endif     
        fprintf(stderr, "usage: %s  [options] 700C|700D|800XA|FSK_LDPC%s InputBinaryDataFile OutputModemRawFile\n"
               "\n"
               "  --testframes N  send N test frames per burst\n"
               "  --bursts     B  send B bursts on N testframes (default 1)\n"
               "  -a amp          maximum amplitude of FSK signal\n"
               "  -c              complex signed 16 bit output format (default real)\n"
               "  --clip  0|1     clipping for reduced PAPR\n"
               "  --txbpf 0|1     bandpass filter\n"
               "  --seq           send packet sequence numbers (breaks testframe BER counting)\n"
               "\n"
               "For FSK_LDPC only:\n\n"
               "  -m      2|4     number of FSK tones\n"        
               "  --Fs    FreqHz  sample rate (default 8000)\n"
               "  --Rs    FreqHz  symbol rate (default 100)\n"
               "  --tone1 FreqHz  freq of first tone (default 1000)\n"
               "  --shift FreqHz  shift between tones (default 200)\n\n"
               , argv[0], f2020);
        fprintf(stderr, "example: $ %s 700D dataBytes.bin samples.s16\n", argv[0]);
        fprintf(stderr, "example: $ %s FSK_LDPC -c --testframes 10 /dev/zero samples.iq16\n\n", argv[0]);
        exit(1);
    }

    use_clip = 0; use_txbpf = 0; testframes = 0; use_complex = 0;
    
    int o = 0;
    int opt_idx = 0;
    while( o != -1 ){
        static struct option long_opts[] = {
            {"testframes", required_argument,  0, 't'},
            {"help",       no_argument,        0, 'h'},
            {"txbpf",      required_argument,  0, 'b'},
            {"clip",       required_argument,  0, 'l'},
            {"Fs",         required_argument,  0, 'f'},
            {"Rs",         required_argument,  0, 'r'},
            {"tone1",      required_argument,  0, '1'},
            {"shift",      required_argument,  0, 's'},
            {"bursts",     required_argument,  0, 'e'},
            {"seq",        no_argument,        0, 'q'},
            {0, 0, 0, 0}
        };
        
        o = getopt_long(argc,argv,"a:ct:hb:l:e:f:r:1:s:m:q",long_opts,&opt_idx);
        
        switch(o) {
        case 'a':
            amp = atof(optarg)/2.0;
            break;
        case 'b':
            use_txbpf = atoi(optarg);
            break;
        case 'c':
            use_complex = 1;
            shorts_per_sample = 2;
            break;
        case 'e':
            Nbursts = atoi(optarg);
            break;
        case 't':
            testframes = 1;
            Nframes = atoi(optarg);
            break;
        case 'l':
            use_clip = atoi(optarg);
            break;
        case 'm':
            adv.M = atoi(optarg);
            break;
        case 'f':
            adv.Fs = atoi(optarg);
            break;
        case 'r':
            adv.Rs = atoi(optarg);
            break;
        case '1':
            adv.first_tone = atoi(optarg);
            break;
        case 's':
            adv.tone_spacing = atoi(optarg);
            break;
        case 'q':
            sequence_numbers = 1;
            break;
        case 'h':
        case '?':
            goto helpmsg;
            break;
        }
    }
    int dx = optind;
    
    if( (argc - dx) < 3) {
        fprintf(stderr, "too few arguments.\n");
        goto helpmsg;
    }

    mode = -1;
    if (!strcmp(argv[dx],"700C")) mode = FREEDV_MODE_700C;
    if (!strcmp(argv[dx],"700D")) mode = FREEDV_MODE_700D;
    #ifdef __LPCNET__
    if (!strcmp(argv[dx],"2020")) mode = FREEDV_MODE_2020;
    #endif
    if (!strcmp(argv[dx],"FSK_LDPC")) mode = FREEDV_MODE_FSK_LDPC;
    if (mode == -1) {
        fprintf(stderr, "Error in mode: %s\n", argv[1]);
        goto helpmsg;        
    }

    if (strcmp(argv[dx+1], "-")  == 0) fin = stdin;
    else if ( (fin = fopen(argv[dx+1],"rb")) == NULL ) {
        fprintf(stderr, "Error opening input file of bytes: %s: %s.\n", argv[dx+1], strerror(errno));
        exit(1);
    }

    if (strcmp(argv[dx+2], "-") == 0) fout = stdout;
    else if ( (fout = fopen(argv[dx+2],"wb")) == NULL ) {
        fprintf(stderr, "Error opening output modem sample file: %s: %s.\n", argv[dx+2], strerror(errno));
        exit(1);
    }

    if (mode != FREEDV_MODE_FSK_LDPC)
        freedv = freedv_open(mode);
    else
        freedv = freedv_open_advanced(mode, &adv);
        
    assert(freedv != NULL);

    /* these are optional ------------------ */
    freedv_set_clip(freedv, use_clip);
    freedv_set_tx_bpf(freedv, use_txbpf);
    freedv_set_tx_amp(freedv, amp);
    
    /* for streaming bytes it's much easier to use modes that have a multiple of 8 payload bits/frame */
    int bytes_per_modem_frame = freedv_get_bits_per_modem_frame(freedv)/8;
    int payload_bytes_per_modem_frame = bytes_per_modem_frame;
    if (mode == FREEDV_MODE_FSK_LDPC) payload_bytes_per_modem_frame -= 2; /* 16 bits used for the CRC */
    fprintf(stderr, "bits_per_modem_frame: %d bytes_per_modem_frame: %d\n", freedv_get_bits_per_modem_frame(freedv), bytes_per_modem_frame);
    assert((freedv_get_bits_per_modem_frame(freedv) % 8) == 0);
    int     n_mod_out = freedv_get_n_nom_modem_samples(freedv);
    uint8_t bytes_in[bytes_per_modem_frame];

    if (mode == FREEDV_MODE_FSK_LDPC) {
        fprintf(stderr, "Frequency: Fs: %4.1f kHz Rs: %4.1f kHz Tone1: %4.1f kHz Shift: %4.1f kHz M: %d \n",
                (float)adv.Fs/1E3, (float)adv.Rs/1E3, (float)adv.first_tone/1E3, (float)adv.tone_spacing/1E3, adv.M);

        if (adv.tone_spacing < adv.Rs) {
            fprintf(stderr, "Need shift: %d > Rs: %d\n", adv.tone_spacing, adv.Rs);
            exit(1);
        }
    }
    
    if ((Nbursts != 1) && (testframes == 0)) {
        fprintf(stderr, "Error: --bursts can only be used with --testframes\n");
        exit(1);
    }
    
    /* optionally set up a known testframe */
    uint8_t testframe_bytes[bytes_per_modem_frame];
    memset(testframe_bytes, 0, bytes_per_modem_frame);
    if (testframes) {
        int bits_per_frame = freedv_get_bits_per_modem_frame(freedv);
        uint8_t testframe_bits[bits_per_frame];
        ofdm_generate_payload_data_bits(testframe_bits, bits_per_frame);
        freedv_pack(testframe_bytes, testframe_bits, bits_per_frame);
    }
    fprintf(stderr, "\n");

    short mod_out_short[2*n_mod_out];
    COMP  mod_out_comp[n_mod_out];
    int frames;

    for(int b=0; b<Nbursts; b++) {

        /* send preamble to help estimators lock up at start of burst */
        int n_preamble = 0;
        if (mode == FREEDV_MODE_FSK_LDPC) {
            if (use_complex == 0) {
                n_preamble = freedv_rawdatapreambletx(freedv, mod_out_short);
            } else {
                n_preamble = freedv_rawdatapreamblecomptx(freedv, mod_out_comp);
                comp_to_short(mod_out_short, mod_out_comp, n_preamble);
            }
            fwrite(mod_out_short, sizeof(short), shorts_per_sample*n_preamble, fout);
        }
        
        /* OK main loop  --------------------------------------- */

        frames = 0;
        while(fread(bytes_in, sizeof(uint8_t), payload_bytes_per_modem_frame, fin) == payload_bytes_per_modem_frame) {
            if (testframes) {
                memcpy(bytes_in, testframe_bytes, bytes_per_modem_frame);
                if (sequence_numbers) bytes_in[0] = (frames+1) & 0xff;
            }
            if (mode == FREEDV_MODE_FSK_LDPC) {
                
                /* This mode requires a CRC in the last two bytes. TODO: consider moving inside freedv_rawdatatx(),
                   although there may be some advantage in leaving the CRC visible to upper layers */
                
                uint16_t crc16 = freedv_gen_crc16(bytes_in, payload_bytes_per_modem_frame);
                bytes_in[bytes_per_modem_frame-2] = crc16 >> 8; 
                bytes_in[bytes_per_modem_frame-1] = crc16 & 0xff; 
            }
            
            if (use_complex == 0) {
                freedv_rawdatatx(freedv, mod_out_short, bytes_in);
            } else {
                freedv_rawdatacomptx(freedv, mod_out_comp, bytes_in);
                comp_to_short(mod_out_short, mod_out_comp, n_mod_out);
            }
            fwrite(mod_out_short, sizeof(short), shorts_per_sample*n_mod_out, fout);
    
            /* if using pipes we don't want the usual buffering to occur */
            if (fout == stdout) fflush(stdout);
            if (fin == stdin) fflush(stdin);

            frames++;       
            if (testframes && (frames >= Nframes)) break;
        }

        /* some silence at the end to allow demod to complete processing */
        
        for(int i=0; i<shorts_per_sample*n_mod_out; i++) mod_out_short[i] = 0;
        fwrite(mod_out_short, sizeof(short), shorts_per_sample*n_mod_out, fout);
        fwrite(mod_out_short, sizeof(short), shorts_per_sample*n_mod_out, fout);
    }
    
    freedv_close(freedv);
    fclose(fin);
    fclose(fout);
    
    return 0;
}

