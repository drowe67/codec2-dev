/*---------------------------------------------------------------------------*\
                                                                             
  FILE........: freedv.c
  AUTHOR......: David Rowe
  DATE CREATED: August 2014
                                                                             
  Functions that implement FreeDV "modes" (or waveforms), useful for
  embedding in other programs.
                                                                             
\*---------------------------------------------------------------------------*/

/*
  Copyright (C) 2012 David Rowe

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

#define FREEDV_1600       0
#define FREEDV_NSAMPLES 320

#include "codec2.h"
#include "codec2_fdmdv.h"
#include "golay23.h"


struct freedv {
    int           mode;
    void         *codec2;
    struct FDMDV *fdmdv
    char         *packed_codec_bits;
    int          *codec_bits;
    int          *tx_bits;
    int          *fdmdv_bits;
    int          *rx_bits;
    int           tx_sync_bit;
};


struct freedv *freedv_open(int mode) {
    struct freedv *f;
    int            Nc, codec2_mode, nbit, nbyte;

    f = (struct freedv*)malloc(sizeof(struct freedv));
    if (f == NULL)
        return NULL;

    f->mode = mode;
    f->tx_sync = 0;

    if (mode == FREEDV_1600) {
        Nc = 16;
        codec2_mode = CODEC2_MODE_1300;
    }

    f->codec2 = codec2_create(Nc);
    if (f->codec2 == NULL)
        return NULL;

    f->fdmdv = fdmdv_create(codec2_mode);
    if (f->fdmdv == NULL)
        return NULL;

    nbit = codec2_bits_per_frame(f->codec2);
    nbyte = (nbit + 7) / 8;
    f->packed_codec_bits = (unsigned char*)malloc(nbyte*sizeof(char));
    f->codec_bits = (int*)malloc(nbit*sizeof(int));

    nbit = 2*fdmdv_bits_per_frame(fdmdv);
    f->tx_bits = (int*)malloc(nbit*sizeof(int));
    f->rx_bits = (int*)malloc(nbit*sizeof(int));
 
    nbit = fdmdv_bits_per_frame(fdmdv);
    f->fdmdv_bits = (int*)malloc(nbit*sizeof(int));

    if ((f->packed_codec_bits == NULL) || (f->codec_bits == NULL) 
        || (f->tx_bits == NULL) || (f->rx_bits == NULL) || (f->fdmdv_bits == NULL))
        return NULL;

    golay23_init();

    return f;
}

void freedv_close(struct freedv *freedv) {
    free(freedv->packed_codec_bits);
    free(freedv->codec_bits);
    free(freedv->tx_bits);
    fdmdv_destroy(freedv->fdmdv);
    codec2_destroy(freedv->codec2);
    free(freedv);
}

int freedv_tx(struct freedv *f, short speech_in[], short mod_out[]) {
    int    bit, byte, i, j;
    int    bits_per_codec_frame, bits_per_fdmdv_frame;
    int    data, codeword1, codeword2;
    COMP   tx_fdm[2*FDMDV_NOM_SAMPLES_PER_FRAME];
    short  tx_fdm_scaled[2*FDMDV_NOM_SAMPLES_PER_FRAME];
     
    bits_per_codec_frame = codec2_bits_per_frame(freedv->codec2);
    bits_per_fdmdv_frame = fdmdv_bits_per_frame(freedv->codec2);

    codec2_encode(f->codec2, f->packet_codec_bits, speech_in);

    /* unpack bits, MSB first */

    bit = 7; byte = 0;
    for(i=0; i<bits_per_input_frame; i++) {
        f->codec_bits[i] = (f->packed_codec_bits[byte] >> bit) & 0x1;
        bit--;
        if (bit < 0) {
            bit = 7;
            byte++;
        }
    }
    
    if (f->mode == MODE_1600) {
            
        /* Protect first 12 out of first 16 excitation bits with (23,12) Golay Code:

           0,1,2,3: v[0]..v[1]
           4,5,6,7: MSB of pitch
           11,12,13,14: MSB of energy

        */

        data = 0;
        for(i=0; i<8; i++) {
            data <<= 1;
            data |= f->codec_bits[i];
        }
        for(i=11; i<15; i++) {
            data <<= 1;
            data |= f->codec_bits[i];
        }
        codeword1 = golay23_encode(data);

        /* now pack output frame with parity bits at end to make them
           as far apart as possible from the data they protect.  Parity
           bits are LSB of the Golay codeword */

        for(i=0; i<bits_per_codec_frame; i++)
            f->tx_bits[i] = f->codec_bits[i];
        for(j=0; i<bits_per_codec_frame+11; i++,j++) {
            f->tx_bits[i] = (codeword1 >> (10-j)) & 0x1;
        }
        f->tx_bits[i] = 0; /* spare bit */
    }

    /* modulate even and odd frames */

    fdmdv_mod(f->fdmdv, tx_fdm, f->tx_bits, &f->tx_sync_bit);
    assert(sync_bit == 1);

    fdmdv_mod(f->fdmdv, &tx_fdm[FDMDV_NOM_SAMPLES_PER_FRAME], &f->tx_bits[bits_per_fdmdv_frame], &f->tx_sync_bit);
    assert(sync_bit == 0);

    for(i=0; i<2*FDMDV_NOM_SAMPLES_PER_FRAME; i++)
        mod_out[i] = FDMDV_SCALE * tx_fdm[i].real;
}

int freedv_nin(struct freedv *f) {
    return f->nin;
}

/* TODO: sync code, SNR threshold (default but can be changed), unit test program */

int freedv_rx(struct freedv *f, short demod_in[], short speech_out[]) {
    COMP  rx_fdm[FDMDV_MAX_SAMPLES_PER_FRAME];
    int   bits_per_codec_frame, bytes_per_codec_frame, bits_per_fdmdv_frame;
    int   reliable_sync_bit, i, j, bit, byte, ;
    int   recd_codeword, codeword1, codeword2;

    bits_per_codec_frame  = codec2_bits_per_frame(freedv->codec2);
    bytes_per_codec_frame = (bits_per_codec_frame + 7) / 8;
    bits_per_fdmdv_frame  = fdmdv_bits_per_frame(freedv->codec2);

    for(i=0; i<f->nin; i++) {
        rx_fdm[i].real = (float)demod_in[i]/FDMDV_SCALE;
        rx_fdm[i].imag = 0;
    }
    nin_prev = nin;
    fdmdv_demod(f->fdmdv, f->fdmdv_bits, &reliable_sync_bit, rx_fdm, &f->nin);
    
    if (reliable_sync_bit == 0) {
        memcpy(f->rx_bits, f->fdmdv_bits, bits_per_fdmdv_frame*sizeof(int));
    }
    else {
        memcpy(&f->rx_bits[bits_per_fdmdv_frame], f->fdmdv_bits, bits_per_fdmdv_frame*sizeof(int));
    }

    if (f->mode == MODE_1600) {
        recd_codeword = 0;
        for(i=0; i<8; i++) {
            recd_codeword <<= 1;
            recd_codeword |= f->rx_bits[i];
        }
        for(i=11; i<15; i++) {
            recd_codeword <<= 1;
            recd_codeword |= f->rx_bits[i];
        }
        for(i=bits_per_codec_frame; i<bits_per_codec_frame+11; i++) {
            recd_codeword <<= 1;
            recd_codeword |= f->rx_bits[i];
        }
        codeword1 = golay23_decode(recd_codeword);
        //codeword1 = recd_codeword;
        //fprintf(stderr, "received codeword1: 0x%x  decoded codeword1: 0x%x\n", recd_codeword, codeword1);
           
        for(i=0; i<bits_per_codec_frame; i++)
            f->codec_bits[i] = f->rx_bits[i];

        for(i=0; i<8; i++) {
            f->codec_bits[i] = (codeword1 >> (22-i)) & 0x1;
        }
        for(i=8,j=11; i<12; i++,j++) {
            f->codec_bits[j] = (codeword1 >> (22-i)) & 0x1;
        }
    }
 
    // pack bits, MSB received first

    bit  = 7;
    byte = 0;
    memset(f->packed_codec_bits, 0,  bytes_per_codec_frame);
    for(i=0; i<bits_per_codec_frame; i++) {
        f->packed_codec_bits[byte] |= (f->codec_bits[i] << bit);
        bit--;
        if(bit < 0) {
            bit = 7;
            byte++;
        }
    }

    codec2_decode(f->codec2, speech_out, f->packed_codec_bits);
}

int freedv_tx_text(struct freedv_mode *mode, char c) {
}

int freedv_rx_text(struct freedv_mode *mode, char *c) {
}

