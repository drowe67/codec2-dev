/*---------------------------------------------------------------------------*\

  FILE........: fsk_mod_ext_vco.c
  AUTHOR......: David Rowe
  DATE CREATED: Feb 2018

  Converts a stream of bits to mFSK raw file of "levels" suitable for
  driving an external VCO, e.g. legacy FM tranmsitter in data mode.
   
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
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

int main(int argc,char *argv[]){
    int   os, m, log2m, i, bit_i, sym;
    float d;
    FILE *fin,*fout;
    
    if(argc<5){
        fprintf(stderr,"usage: %s MbitsPerFSKsymbol OutputSamplesPerSymbol deviationPerlevel InputOneBitPerCharFile OutputVcoRawFile\n",argv[0]);
        exit(1);
    }
    
    /* Extract parameters */
    m = atoi(argv[1]); log2m = log2(m);
    printf("log2m: %d\n", log2m);
    os = atoi(argv[2]);
    d = atof(argv[3]);
    
    uint8_t tx_bits[log2m];
    int16_t rawbuf[os];

    if (strcmp(argv[4],"-")==0){
        fin = stdin;
    } else {
        fin = fopen(argv[4],"r");
    }
	
    if(strcmp(argv[5],"-")==0){
        fout = stdout;
    } else {
        fout = fopen(argv[5],"w");
    }
    
    /* Modulate m bits to levels to drive external VCO */

    while( fread(tx_bits, sizeof(uint8_t), log2m, fin) == log2m ){

        /* generate the symbol number from the bit stream, 
           e.g. 0,1 for 2FSK, 0,1,2,3 for 4FSK */

        sym = bit_i = 0;
        for( i=m; i>>=1; ){
            //fprintf(stderr, "tx_bits[%d] = %d\n", i, tx_bits[bit_i]);
            uint8_t bit = tx_bits[bit_i];
            bit = (bit==1)?1:0;
            sym = (sym<<1)|bit;
            bit_i++;
        }
        //fprintf(stderr, "sym = %d\n", sym);

        /* map 'sym' to VCO drive signal symmetrically about 0,
           separate tones by constant "d"  */
        /* 2 FSK -d/2, +d/2                */
        /* 4 FSK -3*d/2, -d/2, +d/2, 3*d/2 */

        /* note: drive is inverted, a higher tone drives VCO voltage lower */

        float symf = sym;
        float level = d*(((float)m-1)*0.5 - symf);
        assert(level <= 32767.0);
        assert(level >= -32768.0);
        short level_short = (short)level;
        //fprintf(stderr, "level = %f level_short = %d\n\n", level, level_short);
        for(i=0; i<os; i++) {
            rawbuf[i] = level_short;
        }
        fwrite(rawbuf, sizeof(int16_t), os, fout);
        
        if(fin == stdin || fout == stdin){
            fflush(fin);
            fflush(fout);
        }
    }
    
    fclose(fin);
    fclose(fout);

    exit(0);
}
