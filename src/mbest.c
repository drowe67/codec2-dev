/*---------------------------------------------------------------------------*\

  FILE........: mbest.c
  AUTHOR......: David Rowe
  DATE CREATED: Jan 2017

  Multistage vector quantiser search algorithm that keeps multiple
  candidates from each stage.

\*---------------------------------------------------------------------------*/

/*
  Copyright David Rowe 2017

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
#include <string.h>

#include "mbest.h"

struct MBEST *mbest_create(int entries) {
    int           i,j;
    struct MBEST *mbest;

    assert(entries > 0);
    mbest = (struct MBEST *)malloc(sizeof(struct MBEST));
    assert(mbest != NULL);

    mbest->entries = entries;
    mbest->list = (struct MBEST_LIST *)malloc(entries*sizeof(struct MBEST_LIST));
    assert(mbest->list != NULL);

    for(i=0; i<mbest->entries; i++) {
	for(j=0; j<MBEST_STAGES; j++)
	    mbest->list[i].index[j] = 0;
	mbest->list[i].error = 1E32;
    }

    return mbest;
}


void mbest_destroy(struct MBEST *mbest) {
    assert(mbest != NULL);
    free(mbest->list);
    free(mbest);
}


/*---------------------------------------------------------------------------*\

  mbest_insert

  Insert the results of a vector to codebook entry comparison. The
  list is ordered in order or error, so those entries with the
  smallest error will be first on the list.

\*---------------------------------------------------------------------------*/

void mbest_insert(struct MBEST *mbest, int index[], float error) {
    int                i, j, found;
    struct MBEST_LIST *list    = mbest->list;
    int                entries = mbest->entries;

    found = 0;
    for(i=0; i<entries && !found; i++)
	if (error < list[i].error) {
	    found = 1;
	    for(j=entries-1; j>i; j--)
		list[j] = list[j-1];
	    for(j=0; j<MBEST_STAGES; j++)
		list[i].index[j] = index[j];
	    list[i].error = error;
	}
}


#ifdef MBEST_PRINT_OUT
static void mbest_print(char title[], struct MBEST *mbest) {
    int i,j;

    fprintf(stderr, "%s\n", title);
    for(i=0; i<mbest->entries; i++) {
	for(j=0; j<MBEST_STAGES; j++)
	    fprintf(stderr, "  %4d ", mbest->list[i].index[j]);
	fprintf(stderr, " %f\n", mbest->list[i].error);
    }
}
#endif


/*---------------------------------------------------------------------------*\

  mbest_search

  Searches vec[] to a codebbook of vectors, and maintains a list of the mbest
  closest matches.

\*---------------------------------------------------------------------------*/

void mbest_search(
		  const float  *cb,     /* VQ codebook to search         */
		  float         vec[],  /* target vector                 */
		  float         w[],    /* weighting vector              */
		  int           k,      /* dimension of vector           */
		  int           m,      /* number on entries in codebook */
		  struct MBEST *mbest,  /* list of closest matches       */
		  int           index[] /* indexes that lead us here     */
)
{
   float   e;
   int     i,j;
   float   diff;

   for(j=0; j<m; j++) {
	e = 0.0;
	for(i=0; i<k; i++) {
	    diff = cb[j*k+i]-vec[i];
	    e += powf(diff*w[i],2.0);
	}
	index[0] = j;
	mbest_insert(mbest, index, e);
   }
}

