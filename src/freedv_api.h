/*---------------------------------------------------------------------------*\
                                                                             
  FILE........: freedv_api.h
  AUTHOR......: David Rowe
  DATE CREATED: August 2014
                                                                             
  Library of API functions that implement FreeDV "modes", useful for
  embedding FreeDV in other programs.
                                                                             
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

#ifndef __FREEDV__


#define FREEDV_MODE_1600    0
#define FREEDV_NSAMPLES   320


struct freedv {
    int            mode;
    void          *codec2;
    struct FDMDV  *fdmdv;
    unsigned char *packed_codec_bits;
    int           *codec_bits;
    int           *tx_bits;
    int           *fdmdv_bits;
    int           *rx_bits;
    int            tx_sync_bit;
    float          snr_thresh;
    int            nin;
};


struct freedv *freedv_open(int mode);
void freedv_close(struct freedv *freedv);
void freedv_tx(struct freedv *f, short mod_out[], short speech_in[]);
int freedv_nin(struct freedv *f);
int freedv_rx(struct freedv *f, short speech_out[], short demod_in[]);


#endif
