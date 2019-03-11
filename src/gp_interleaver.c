/*---------------------------------------------------------------------------*\

  FILE........: gp_interleaver.c
  AUTHOR......: David Rowe
  DATE CREATED: April 2018

  Golden Prime Interleaver. My interpretation of "On the Analysis and
  Design of Good Algebraic Interleavers", Xie et al,eq (5).

  See also octave/gp_interleaver.m

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
#include <stdio.h>
#include "gp_interleaver.h"

/*
  Choose b for Golden Prime Interleaver.  b is chosen to be the
  closest integer, which is relatively prime to N, to the Golden
  section of N.

  Implemented with a LUT in C for convenience, Octave version
  has a more complete implementation.  If you find you need some more
  numbers head back to the Octave choose_interleaver_b() function.
*/

int b_table[] = {
  112,71,
  210,131,
  224,139,
  252,157,
  420,263,
  448,277,
  504,313,
  672,419,
  896,557,
  1120,701,
  1344,839,
  1568,971,
  1792,1109,
  2016,1249,
  2240,1399,
  2464,1523,
  2688,1663,
  2912,1801,
  3136,1949,
  3360,2081,
  3584,2213
};

int choose_interleaver_b(int Nbits)
{
    int i;

    for(i=0; i<sizeof(b_table)/(2*sizeof(int)); i+=2) {
        if (b_table[i] == Nbits) {
            return b_table[i+1];
        }
    }

    /* if we get to here it means a Nbits we don't have in our table so choke */
    
    fprintf(stderr, "Nbits: %d, b not found!\n", Nbits);
    assert(0);
}


void gp_interleave_comp(COMP interleaved_frame[], COMP frame[], int Nbits) {
  int b = choose_interleaver_b(Nbits);
  int i,j;
  for (i=0; i<Nbits; i++) {
    j = (b*i) % Nbits;
    interleaved_frame[j] = frame[i];
  }
}

void gp_deinterleave_comp(COMP frame[], COMP interleaved_frame[], int Nbits) {
  int b = choose_interleaver_b(Nbits);
  int i,j;
  for (i=0; i<Nbits; i++) {
    j = (b*i) % Nbits;
    frame[i] =  interleaved_frame[j];
  }
}

void gp_interleave_float(float interleaved_frame[], float frame[], int Nbits) {
  int b = choose_interleaver_b(Nbits);
  int i,j;
   
  for (i=0; i<Nbits; i++) {
    j = (b*i) % Nbits;
    interleaved_frame[j] = frame[i];
  }
}

void gp_deinterleave_float(float frame[], float interleaved_frame[], int Nbits) {
  int b = choose_interleaver_b(Nbits);
  int i,j;
   
  for (i=0; i<Nbits; i++) {
    j = (b*i) % Nbits;
    frame[i] = interleaved_frame[j];
  }
}

