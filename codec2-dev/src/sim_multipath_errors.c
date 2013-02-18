/*---------------------------------------------------------------------------*\

  FILE........: sim_multipath_errors.c
  AUTHOR......: David Rowe
  DATE CREATED: 13/2/2013

  Simulates error bursts on each FDMDV carrier from HF multipath.  Useful for
  testing error protection schemes.

\*---------------------------------------------------------------------------*/

/*
  Copyright (C) 2010 David Rowe

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

#include "codec2.h"

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>

#define SYMBOL_PERIOD 0.02
   
int iterate_state(int    state,
                  int    carrier,
                  float  t, 
                  float  burst_length, 
                  float  burst_offset)
{
    int next_state = state;

    //printf("carrier: %2d state: %d t: %4.3f offset: %f\n", carrier, state, t, burst_offset);

    switch(state) {
    case 0:

        /* clear channel state - no bit errors */

        /* start of error burst in each carrier is shifted by
           burst_offset wrt previous carrier */

        if (t > burst_offset*(carrier+1))
            next_state = 1;
        break;

    case 1:
                
        /* burst error state - 50% bit error rate */

        if (t > (burst_offset*(carrier+1) + burst_length))
            next_state = 2;
        break;

    }
 
    return next_state;
}


int main(int argc, char *argv[])
{
    FILE          *fin;
    FILE          *fout;
    unsigned short bits_in, bit;
    int            i, nbits_in, nbits_out, n_carriers, bits_per_carrier;
    unsigned char  byte, bits_out;
    int           *state, carrier;
    float          t, burst_length, burst_period, burst_offset, r;
    int            errors,bits;

    if (argc < 7) {
	printf("%s InputBitFile OutputBitFile nCarriers BitsPerCarrier "
               "burstPeriod burstLength nextCarrierOffset\n", argv[0]);
	exit(1);
    }

    if (strcmp(argv[1], "-")  == 0) fin = stdin;
    else if ( (fin = fopen(argv[1],"rb")) == NULL ) {
	fprintf(stderr, "Error opening input bit file: %s: %s.\n",
         argv[1], strerror(errno));
	exit(1);
    }

    if (strcmp(argv[2], "-") == 0) fout = stdout;
    else if ( (fout = fopen(argv[2],"wb")) == NULL ) {
	fprintf(stderr, "Error opening output speech file: %s: %s.\n",
         argv[2], strerror(errno));
	exit(1);
    }

    n_carriers = atoi(argv[3]);
    bits_per_carrier = atoi(argv[4]);
    
    burst_period = atof(argv[5]);
    burst_length = atof(argv[6]);
    burst_offset = atof(argv[7]);

    state = (int*)malloc(sizeof(int)*n_carriers);
    assert(state != NULL);
    for(i=0; i<n_carriers; i++)
        state[i] = 0;

    /* We have packed bits for input and output, which must be read as
       bytes.  We would like to process these packed bits at the modem
       frame rate.  e.g. for 1400 bit/s, the modem frame rate is 28
       bit/frame of 3.5 bytes.

       while !feof
         Read a byte
         while enough bits left
           apply burst error
           
    */

    t = 0.0;
    carrier = 0;
    nbits_in = 0;
    bits_in = 0;
    nbits_out = 0;
    bits_out = 0;
    errors = bits = 0;
    while(fread(&byte, sizeof(char), 1, fin) == 1) {

        /* insert latest byte into bits buffer, MSB is first in time order */
        
        bits_in <<= 8;
        bits_in |= (short)byte;
        nbits_in += 8;

        //printf("byte: 0x%x bits_in: 0x%x nbits_in: %d\n", byte, bits_in, nbits_in);

        while(nbits_in >= bits_per_carrier) {

            /* iterate error state for every carrier to see if this carrier is 
               knocked out by a simulated fade */

            state[carrier] = iterate_state(state[carrier], 
                                           carrier,
                                           t, 
                                           burst_length, burst_offset);

            /* apply error model to bits in this carrier */

            for(i=0; i<bits_per_carrier; i++) {

                bit = (bits_in >> (nbits_in - 1)) & 0x1;               
                nbits_in--;
                //printf("bit: %d\n", bit);

                bits++;
                if (state[carrier] == 1) {
                    r = (float)rand()/RAND_MAX;
                    if (r < 0.5) {
                        bit ^= 1;
                        errors++;
                        fprintf(stderr, "e");
                    }
                    else
                        fprintf(stderr, ".");
                }
                else
                    fprintf(stderr, ".");

                /* pack processed bits and fwrite when we get a byte */

                bits_out |= bit;
                nbits_out++;
                if (nbits_out == 8) {
                    //printf("bits_out: 0x%x\n", bits_out);
                    fwrite(&bits_out, sizeof(char), 1, fout); 
                    nbits_out = 0;
                    bits_out = 0;
                }
                else
                    bits_out <<= 1;

            }

            /* move to next carrier */

            carrier++;
            fprintf(stderr, " ");
            if (carrier == n_carriers) {
                carrier = 0;
                t += SYMBOL_PERIOD;
                if (t > burst_period) {
                    t = 0.0;
                    for(i=0; i<n_carriers; i++)
                        state[i] = 0;
                }
                fprintf(stderr, "\n");
            }
        }
    }

    fclose(fin);
    fclose(fout);
    free(state);

    fprintf(stderr,"ber: %4.3f\n", (float)errors/bits);

    return 0;
}
