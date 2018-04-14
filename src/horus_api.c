/*---------------------------------------------------------------------------*\

  FILE........: horus_api.c
  AUTHOR......: David Rowe
  DATE CREATED: March 2018

  Library of API functions that implement High Altitude Balloon (HAB)
  telemetry modems and protocols for Project Horus.  May also be useful for
  other HAB projects.

\*---------------------------------------------------------------------------*/

/*
  Copyright (C) 2018 David Rowe

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

#include "horus_api.h"
#include "fsk.h"
#include "horus_l2.h"

#define MAX_UW_LENGTH                  100
#define HORUS_API_VERSION                1    /* unique number that is bumped if API changes */
#define HORUS_BINARY_NUM_BITS          360    /* fixed number of bytes in binary payload     */
#define HORUS_BINARY_NUM_PAYLOAD_BYTES  22    /* fixed number of bytes in binary payload     */

struct horus {
    int         mode;
    int         verbose;
    struct FSK *fsk;                 /* states for FSK modem                */
    int         Fs;                  /* sample rate in Hz                   */
    int         mFSK;                /* number of FSK tones                 */
    int         Rs;                  /* symbol rate in Hz                   */
    int8_t      uw[MAX_UW_LENGTH];   /* unique word bits mapped to +/-1     */
    int         uw_thresh;           /* threshold for UW detection          */
    int         uw_len;              /* length of unique word               */
    int         max_packet_len;      /* max length of a telemetry packet    */
    uint8_t    *rx_bits;             /* buffer of received bits             */
    int         rx_bits_len;         /* length of rx_bits buffer            */
    int         crc_ok;              /* most recent packet checksum results */
};

/* Unique word for Horus RTTY 7 bit '$' character, 3 sync bits,
   repeated 5 times */

int8_t uw_horus_rtty[] = {
  0,0,1,0,0,1,0,1,1,0,
  0,0,1,0,0,1,0,1,1,0,
  0,0,1,0,0,1,0,1,1,0,
  0,0,1,0,0,1,0,1,1,0,
  0,0,1,0,0,1,0,1,1,0
};

/* Unique word for Horus Binary */

int8_t uw_horus_binary[] = {
    0,0,1,0,0,1,0,0,
    0,0,1,0,0,1,0,0    
};

struct horus *horus_open (int mode) {
    int i;
    assert((mode == HORUS_MODE_RTTY) || (mode == HORUS_MODE_BINARY));

    struct horus *hstates = (struct horus *)malloc(sizeof(struct horus));
    assert(hstates != NULL);

    hstates->Fs = 48000; hstates->Rs = 100; hstates->verbose = 0; hstates->mode = mode;

    if (mode == HORUS_MODE_RTTY) {
        hstates->mFSK = 2;
        hstates->max_packet_len = 1000;

        /* map UW to make it easier to search for */

        for (i=0; i<sizeof(uw_horus_rtty); i++) {
            hstates->uw[i] = 2*uw_horus_rtty[i] - 1;
        }        
        hstates->uw_len = sizeof(uw_horus_rtty);
        hstates->uw_thresh = sizeof(uw_horus_rtty) - 2;  /* allow a few bit errors in UW detection */
        hstates->rx_bits_len = hstates->max_packet_len;
    }

    if (mode == HORUS_MODE_BINARY) {
        hstates->mFSK = 4;
        hstates->max_packet_len = HORUS_BINARY_NUM_BITS;
        for (i=0; i<sizeof(uw_horus_binary); i++) {
            hstates->uw[i] = 2*uw_horus_binary[i] - 1;
        }
        hstates->uw_len = sizeof(uw_horus_binary);
        hstates->uw_thresh = sizeof(uw_horus_binary) - 2; /* allow a few bit errors in UW detection */
        horus_l2_init();
        hstates->rx_bits_len = hstates->max_packet_len;
    }
   
    hstates->fsk = fsk_create(hstates->Fs, hstates->Rs, hstates->mFSK, 1000, 2*hstates->Rs);

    /* allocate enough room for two packets so we know there will be
       one complete packet if we find a UW at start */
    
    hstates->rx_bits_len += hstates->fsk->Nbits;
    hstates->rx_bits = (uint8_t*)malloc(hstates->rx_bits_len);
    assert(hstates->rx_bits != NULL);
    for(i=0; i<hstates->rx_bits_len; i++) {
        hstates->rx_bits[i] = 0;
    }

    hstates->crc_ok = 0;
    
    return hstates;
}

void horus_close (struct horus *hstates) {
    assert(hstates != NULL);
    fsk_destroy(hstates->fsk);
    free(hstates->rx_bits);
    free(hstates);
}

uint32_t horus_nin(struct horus *hstates) {
    assert(hstates != NULL);
    int nin = fsk_nin(hstates->fsk);
    assert(nin <= horus_get_max_demod_in(hstates));
    return nin;
}

int horus_find_uw(struct horus *hstates, int n) {
    int i, j, corr, mx, mx_ind;
    int rx_bits_mapped[n+hstates->uw_len];
    
    /* map rx_bits to +/-1 for UW search */

    for(i=0; i<n+hstates->uw_len; i++) {
        rx_bits_mapped[i] = 2*hstates->rx_bits[i] - 1;
    }
    
    /* look for UW  */

    mx = 0; mx_ind = 0;
    for(i=0; i<n; i++) {

        /* calculate correlation between bit stream and UW */
        
        corr = 0;
        for(j=0; j<hstates->uw_len; j++) {
            corr += rx_bits_mapped[i+j]*hstates->uw[j];
        }

        /* peak pick maximum */
        
        if (corr > mx) {
            mx = corr;
            mx_ind = i;
        }
    }

    if (mx >= hstates->uw_thresh) {
        return mx_ind;
    } else {
        return -1;
    }
}

int hex2int(char ch) {
    if (ch >= '0' && ch <= '9')
        return ch - '0';
    if (ch >= 'A' && ch <= 'F')
        return ch - 'A' + 10;
    if (ch >= 'a' && ch <= 'f')
        return ch - 'a' + 10;
    return -1;
}


int extract_horus_rtty(struct horus *hstates, char ascii_out[], int uw_loc) {
    const int nfield = 7;                               /* 7 bit ASCII                    */
    const int npad   = 3;                               /* 3 sync bits between characters */
    int st = uw_loc;                                    /* first bit of first char        */
    int en = hstates->max_packet_len - nfield;          /* last bit of max length packet  */

    int      i, j, endpacket, nout, crc_ok;
    uint8_t  char_dec;
    char    *pout, *ptx_crc;
    uint16_t rx_crc, tx_crc;

    pout = ascii_out; nout = 0; crc_ok = 0; endpacket = 0; rx_crc = tx_crc = 0;
    
    for (i=st; i<en; i+=nfield+npad) {

        /* assemble char LSB to MSB */

        char_dec = 0;
        for(j=0; j<nfield; j++) {
            assert(hstates->rx_bits[i+j] <= 1);
            char_dec |= hstates->rx_bits[i+j] * (1<<j);
        }
        if (hstates->verbose) {
            fprintf(stderr, "i: %4d 0x%02x %c ", i, char_dec, char_dec);
            if ((nout % 6) == 0) {
                fprintf(stderr, "\n");
            }
        }

        /*  if we find a '*' that's the end of the packet for RX CRC calculations */

        if (!endpacket && (char_dec == 42)) {
            endpacket = 1;
            rx_crc = horus_l2_gen_crc16((uint8_t*)&ascii_out[5], nout-5);
            ptx_crc = pout + 1; /* start of tx CRC */
        }

        /* build up output array, really only need up to tx crc but
           may end up going further */
        
        *pout++ = (char)char_dec;
        nout++;
        
    }

    /* if we found the end of packet flag and have enough chars to compute checksum ... */

    //fprintf(stderr, "\n\ntx CRC...\n");
    if (endpacket && (pout > (ptx_crc+3))) {
        tx_crc = 0;
        for(i=0; i<4; i++) {
            tx_crc <<= 4;
            tx_crc |= hex2int(ptx_crc[i]);
            //fprintf(stderr, "ptx_crc[%d] %c 0x%02X tx_crc: 0x%04X\n", i, ptx_crc[i], hex2int(ptx_crc[i]), tx_crc);
        }
        crc_ok = (tx_crc == rx_crc);
        *(ptx_crc+4) = 0;  /* terminate ASCII string */
    }
    else {
        *ascii_out = 0;
    }

    if (hstates->verbose) {
        fprintf(stderr, "\n endpacket: %d nout: %d tx_crc: 0x%04x rx_crc: 0x%04x\n",
                endpacket, nout, tx_crc, rx_crc);
    }
    
    /* make sure we don't overrun storage */
    
    assert(nout <= horus_get_max_ascii_out_len(hstates));

    hstates->crc_ok = crc_ok;
    
    return crc_ok;
}


int extract_horus_binary(struct horus *hstates, char hex_out[], int uw_loc) {
    const int nfield = 8;                      /* 8 bit binary                   */
    int st = uw_loc;                           /* first bit of first char        */
    int en = uw_loc + hstates->max_packet_len; /* last bit of max length packet  */

    int      j, b, nout;
    uint8_t  rxpacket[hstates->max_packet_len];
    uint8_t  rxbyte, *pout;
 
    /* convert bits to a packet of bytes */
    
    pout = rxpacket; nout = 0;
    
    for (b=st; b<en; b+=nfield) {

        /* assemble bytes MSB to LSB */

        rxbyte = 0;
        for(j=0; j<nfield; j++) {
            assert(hstates->rx_bits[b+j] <= 1);
            rxbyte <<= 1;
            rxbyte |= hstates->rx_bits[b+j];
        }
        
        /* build up output array */
        
        *pout++ = rxbyte;
        nout++;
    }

    if (hstates->verbose) {
        fprintf(stderr, "nout: %d\nReceived Packet before decoding:\n", nout);
        for (b=0; b<nout; b++) {
            fprintf(stderr, "%02X", rxpacket[b]);
        }
        fprintf(stderr, "\n");
    }
    
    uint8_t payload_bytes[HORUS_BINARY_NUM_PAYLOAD_BYTES];
    horus_l2_decode_rx_packet(payload_bytes, rxpacket, HORUS_BINARY_NUM_PAYLOAD_BYTES);

    uint16_t crc_tx, crc_rx;
    crc_rx = horus_l2_gen_crc16(payload_bytes, HORUS_BINARY_NUM_PAYLOAD_BYTES-2);
    crc_tx = *(uint16_t*)&payload_bytes[HORUS_BINARY_NUM_PAYLOAD_BYTES-2];
    
    if (hstates->verbose) {
        fprintf(stderr, "crc_tx: %04X crc_rx: %04X\n", crc_tx, crc_rx);
    }
    
    /* convert to ASCII string of hex characters */

    hex_out[0] = 0;
    char hex[3];
    for (b=0; b<HORUS_BINARY_NUM_PAYLOAD_BYTES; b++) {
        sprintf(hex, "%02X", payload_bytes[b]);
        strcat(hex_out, hex);
    }
   
    if (hstates->verbose) {
        fprintf(stderr, "nout: %d\nDecoded Payload bytes:\n%s", nout, hex_out);
    }

    hstates->crc_ok = (crc_tx == crc_rx);

    /* binary packets always marked as OK, as next layer determines validity */
    
    return 1;
}


int horus_rx(struct horus *hstates, char ascii_out[], short demod_in[]) {
    int i, j, uw_loc, packet_detected;
    
    assert(hstates != NULL);
    packet_detected = 0;

    int Nbits = hstates->fsk->Nbits;
    int rx_bits_len = hstates->rx_bits_len;
    
    if (hstates->verbose) {
        fprintf(stderr, "max_packet_len: %d rx_bits_len: %d Nbits: %d\n",
                hstates->max_packet_len, rx_bits_len, Nbits);
    }
    
    /* shift buffer of bits to make room for new bits */

    for(i=0,j=Nbits; j<rx_bits_len; i++,j++) {
        hstates->rx_bits[i] = hstates->rx_bits[j];
    }
                   
    /* demodulate latest bits */

    COMP demod_in_comp[hstates->fsk->nin];
    for (i=0; i<hstates->fsk->nin; i++) {
        demod_in_comp[i].real = demod_in[i];
        demod_in_comp[i].imag = 0;
    }
    fsk_demod(hstates->fsk, &hstates->rx_bits[rx_bits_len-Nbits], demod_in_comp);
    
    /* UW search to see if we can find the start of a packet in the buffer */
    
    if ((uw_loc = horus_find_uw(hstates, Nbits)) != -1) {

        if (hstates->verbose) {
            fprintf(stderr, "uw_loc: %d mode: %d\n", uw_loc, hstates->mode);
        }
        
        /* OK we have found a unique word, and therefore the start of
           a packet, so lets try to extract valid packets */

        if (hstates->mode == HORUS_MODE_RTTY) {
            packet_detected = extract_horus_rtty(hstates, ascii_out, uw_loc);
        }
        if (hstates->mode == HORUS_MODE_BINARY) {
            packet_detected = extract_horus_binary(hstates, ascii_out, uw_loc);
            //#define DUMP_BINARY_PACKET
            #ifdef DUMP_BINARY_PACKET
            FILE *f = fopen("packetbits.txt", "wt"); assert(f != NULL);
            for(i=0; i<hstates->max_packet_len; i++) {
                fprintf(f,"%d ", hstates->rx_bits[uw_loc+i]);
            }
            fclose(f);
            exit(0);
            #endif
        }
    }
     
    return packet_detected;
}

int horus_get_version(void) {
    return HORUS_API_VERSION;
}

int horus_get_mode(struct horus *hstates) {
    assert(hstates != NULL);
    return hstates->mode;
}

int horus_get_Fs(struct horus *hstates) {
    assert(hstates != NULL);
    return hstates->Fs;
}

int horus_get_mFSK(struct horus *hstates) {
    assert(hstates != NULL);
    return hstates->mFSK;
}

int horus_get_max_demod_in(struct horus *hstates) {
    /* copied from fsk_demod.c, a nicer fsk_max_nin function would be useful */
    return sizeof(short)*(hstates->fsk->N + hstates->fsk->Ts*2);
}

int horus_get_max_ascii_out_len(struct horus *hstates) {
    assert(hstates != NULL);
    if (hstates->mode == HORUS_MODE_RTTY) {
        return hstates->max_packet_len/10;     /* 7 bit ASCII, plus 3 sync bits */
    }
    if (hstates->mode == HORUS_MODE_BINARY) {
        return HORUS_BINARY_NUM_PAYLOAD_BYTES;
    }
    assert(0); /* should never get here */
    return 0;
}

void horus_get_modem_stats(struct horus *hstates, int *sync, float *snr_est) {
    struct MODEM_STATS stats;
    assert(hstates != NULL);

    /* TODO set sync if UW found "recently", but WTF is recently? Maybe need a little state 
       machine to "blink" sync when we get a packet */

    *sync = 0;
    
    /* SNR scaled from Eb/No est returned by FSK to SNR in 3000 Hz */

    fsk_get_demod_stats(hstates->fsk, &stats);
    *snr_est = stats.snr_est + 10*log10(hstates->Rs/3000);
}

void horus_get_modem_extended_stats (struct horus *hstates, struct MODEM_STATS *stats) {
    int i;
    
    assert(hstates != NULL);

    fsk_get_demod_stats(hstates->fsk, stats);

    stats->snr_est = stats->snr_est + 10*log10(hstates->Rs/3000);

    assert(hstates->mFSK <= MODEM_STATS_MAX_F_EST);
    for (i=0; i<hstates->mFSK; i++) {
        stats->f_est[i] = hstates->fsk->f_est[i];
    }
}

void horus_set_verbose(struct horus *hstates, int verbose) {
    assert(hstates != NULL);
    hstates->verbose = verbose;
}

int horus_crc_ok(struct horus *hstates) {
    assert(hstates != NULL);
    return hstates->crc_ok;
}

