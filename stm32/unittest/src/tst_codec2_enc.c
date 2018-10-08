/*---------------------------------------------------------------------------*\

  FILE........: tst_codec2_enc.c, (derived from codec2_profile.c)
  AUTHOR......: David Rowe, Don Reid
  DATE CREATED: 30 May 2013, Oct 2018

  Test Codec 2 operation on the STM32F4.

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

/* This is a unit test implementation of the Codec2_encode function.
 *
 * Typical run:

    # Copy a raw audio file to stm_in.raw. 
    cp ../../../../raw/hts1.raw stm_in.raw

    # Run x86 command for reference output
    c2enc 1600 stm_in.raw mod_ref.raw 

    # Create config
    echo "20000000" > stm_cfg.txt

    # Run stm32
    run_stm32_prog ../../src/tst_codec2_enc.elf --load

    # Compare outputs
    comare_ints -b 1 mod_ref.raw stm_out.raw

    # Manual play (and listen)
    c2dec 1600 mod_ref.raw ref_out.raw
    aplay -f S16_LE  ref_out.raw
    #
    c2dec 1600 stm_out.raw stm_dec.raw
    aplay -f S16_LE stm_dec.raw

 */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>

#include "codec2.h"

#include "stm32f4xx_conf.h"
#include "stm32f4xx.h"
#include "semihosting.h"
#include "machdep.h"

int main(int argc, char *argv[]) {
    FILE           *fcfg, *fin, *fout;

    struct CODEC2  *codec2;
    short          *inbuf, *outbuf;
    unsigned char  *bits;
    int             nsam, nbit, nbyte;
    int             frame;

    ////////
    // Semihosting
    semihosting_init();

    ////////
    // Test configuration, read from stm_cfg.txt
    int     config_mode;        // 0
    //int     config_verbose;     // 6
    //int     config_profile;     // 7
    char config[8];
    fcfg = fopen("stm_cfg.txt", "r");
    if (fcfg == NULL) {
        fprintf(stderr, "Error opening config file\n");
        exit(1);
    }
    if (fread(&config[0], 1, 8, fcfg) != 8) {
        fprintf(stderr, "Error reading config file\n");
        exit(1);
    }
    config_mode = config[0] - '0';
    //config_verbose = config[6] - '0';
    //config_profile = config[7] - '0';
    fclose(fcfg);

    ////////
    // Setup
    codec2 = codec2_create(config_mode);
    nsam = codec2_samples_per_frame(codec2);
    outbuf = (short*)malloc(nsam*sizeof(short));
    inbuf = (short*)malloc(nsam*sizeof(short));
    nbit = codec2_bits_per_frame(codec2);
    nbyte = (nbit + 7) / 8;
    bits = (unsigned char*)malloc(nbyte*sizeof(char));

    fin = fopen("stm_in.raw", "rb");
    if (fin == NULL) {
        printf("Error opening input file\n");
        exit(1);
    }

    fout = fopen("stm_out.raw", "wb");
    if (fout == NULL) {
        printf("Error opening output file\n");
        exit(1);
    }
   
    #ifdef DUMP
    dump_on("stm32f4");
    #endif
    frame = 0;

    while (fread(inbuf, sizeof(short), nsam, fin) == nsam) {
        //PROFILE_SAMPLE(enc_start);
        codec2_encode(codec2, bits, inbuf);
        //PROFILE_SAMPLE_AND_LOG2(, enc_start, "  enc");
        fwrite(bits, sizeof(char), nbyte, fout);
        printf("frame: %d\n", ++frame);
        //machdep_profile_print_logged_samples();
    }

    #ifdef DUMP
    dump_off("sm32f4");
    #endif

    fclose(fin);
    fclose(fout);
    free(inbuf);
    free(outbuf);
    free(bits);
    codec2_destroy(codec2);

    printf("\nEnd of Test\n");
    fclose(stdout);
    fclose(stderr);

    return(0);
}

/* vi:set ts=4 et sts=4: */
