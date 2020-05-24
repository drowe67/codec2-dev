/*---------------------------------------------------------------------------*\

  FILE........: tfreedv_mode_800XA.c
  AUTHOR......: Jeroen Vreeken
  DATE CREATED: 24 May 2020

  Tests specific for mode 800XA

\*---------------------------------------------------------------------------*/
/*
  Copyright (C) 2020 Jeroen Vreeken <jeroen@vreeken.net>

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 2.1 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.


  You should have received a copy of the GNU Lesser General Public License
  along with this program; if not, see <http://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include "freedv_api.h"

int main(int argc, char **argv)
{
    struct freedv *f;
    int i;

    printf("freedv_api tests for mode 800XA\n");

    printf("freedv_open(FREEDV_MODE_800XA) ");
    f = freedv_open(FREEDV_MODE_800XA);
    if (!f) {
        printf("Failed to open\n");
        goto fail;
    }
    printf("Passed\n");

    printf("freedv_get_mode() ");
    int mode = freedv_get_mode(f);
    if (mode != FREEDV_MODE_800XA) {
        printf("mode %d does not match FREEDV_MODE_800XA %d\n", mode, FREEDV_MODE_800XA);
        goto fail;
    }
    printf("Passed\n");

    printf("freedv_get_n_max_modem_samples() ");
    int max_samples = freedv_get_n_max_modem_samples(f);
    if (max_samples != 660) {
        printf("modem max samples %d != 660\n", max_samples);
        goto fail;
    }
    printf("%d Passed\n", max_samples);

    printf("freedv_get_n_nom_modem_samples() ");
    int nom_samples = freedv_get_n_nom_modem_samples(f);
    if (nom_samples != 640) {
        printf("modem nom samples %d != 640\n", nom_samples);
        goto fail;
    }
    printf("%d Passed\n", nom_samples);

    printf("freedv_get_n_speech_samples() ");
    int speech_samples = freedv_get_n_speech_samples(f);
    if (speech_samples != 640) {
        printf("Expected 640 speech samples, got %d\n", speech_samples);
        goto fail;
    }
    printf("%d Passed\n", speech_samples);

    printf("freedv_get_n_bits_per_codec_frame() ");
    int codec_bits = freedv_get_bits_per_codec_frame(f);
    if (codec_bits != 28) {
        printf("Expected 28 codec bits, got %d\n", codec_bits);
	goto fail;
    }
    printf("%d Passed\n", codec_bits);

    printf("freedv_get_n_bits_per_modem_frame() ");
    int frame_bits = freedv_get_bits_per_modem_frame(f);
    if (frame_bits != 56) {
        printf("Expected 56 codec bits, got %d\n", frame_bits);
	goto fail;
    }
    printf("%d Passed\n", frame_bits);

    printf("freedv_rawdatatx()/freedv_rawdatarx() ");
    int frames = 0;
    int fails = 0;
    {
        short mod[nom_samples * 10];
	/* Note: A codec frame is only 3.5 bytes!
	   so th fourth and eight bytes will be half empty!
	 */
        unsigned char payload[8] = { 0x11, 0x22, 0x33, 0x40, 0x55, 0x66, 0x77, 0x80 };
        for (i = 0; i < 10; i ++) {
	    freedv_rawdatatx(f, mod + i * nom_samples, payload);
        }
        int nin = 0;
        for (i = 0; i < nom_samples * 9; i += nin) {
            nin = freedv_nin(f);
	    unsigned char rx_payload[8] = {0};
            int r = freedv_rawdatarx(f, rx_payload, mod + i);
            if (r) {
	        int b;
                for (b = 0; b < 8; b++) {
	    	    if (payload[b] != rx_payload[b]) {
		        printf("Received codec bits 0x%02x do not match expected 0x%02x\n", rx_payload[b], payload[b]);
		        fails++;
                    }
                }
	        frames++;
	    }
        }
    }
    if (!frames) {
    	printf("Did not decode any frames successfully\n");
	goto fail;
    }

    printf("Tests passed\n");
    return 0;
fail:
    printf("Test failed\n");
    return 1;
}
