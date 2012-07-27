/*--------------------------------------------------------------------------*\

	FILE........: vqtrainph.c
	AUTHOR......: David Rowe
	DATE CREATED: 27 July 2012

	This program trains phase vector quantisers.  Modified from
	vqtrain.c

\*--------------------------------------------------------------------------*/

/*
  Copyright (C) 2012 David Rowe

  All rights reserved.

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU Lesser General Public License version 2, as
  published by the Free Software Foundation.  This program is
  distributed in the hope that it will be useful, but WITHOUT ANY
  WARRANTY; without even the implied warranty of MERCHANTABILITY or
  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public
  License for more details.

  You should have received a copy of the GNU Lesser General Public License
  along with this program; if not, see <http://www.gnu.org/licenses/>.
*/

/*-----------------------------------------------------------------------*\

				INCLUDES

\*-----------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <assert.h>

typedef struct {
    float real;
    float imag;
} COMP;

/*-----------------------------------------------------------------------* \

				DEFINES

\*-----------------------------------------------------------------------*/

#define	DELTAQ 	0.01		/* quiting distortion			*/
#define	MAX_STR	80		/* maximum string length		*/

/*-----------------------------------------------------------------------*\

			FUNCTION PROTOTYPES

\*-----------------------------------------------------------------------*/

void zero(COMP v[], int k);
void acc(COMP v1[], COMP v2[], int k);
void norm(COMP v[], int k, int n);
int quantise(COMP cb[], COMP vec[], int k, int m, float *se);
void print_vec(COMP cb[], int k, int m);

/*-----------------------------------------------------------------------* \

				MAIN

\*-----------------------------------------------------------------------*/

int main(int argc, char *argv[]) {
    int    k,m;		/* dimension and codebook size			*/
    COMP   *vec;	/* current vector 				*/
    COMP   *cb;		/* vector codebook				*/
    COMP   *cent;	/* centroids for each codebook entry		*/
    int    *n;		/* number of vectors in this interval		*/
    int     J;		/* number of vectors in training set		*/
    int     ind;	/* index of current vector			*/
    float   se;		/* squared error for this iteration		*/
    float   Dn,Dn_1;	/* current and previous iterations distortion	*/
    float   delta;	/* improvement in distortion 			*/
    FILE   *ftrain;	/* file containing training set			*/
    FILE   *fvq;	/* file containing vector quantiser		*/
    int     ret;
    int     i,j;

    /* Interpret command line arguments */

    if (argc != 5)	{
	printf("usage: %s TrainFile K(dimension) M(codebook size) VQFile\n", argv[0]);
	exit(1);
    }

    /* Open training file */

    ftrain = fopen(argv[1],"rb");
    if (ftrain == NULL) {
	printf("Error opening training database file: %s\n",argv[1]);
	exit(1);
    }

    /* determine k and m, and allocate arrays */

    k = atol(argv[2]);
    m = atol(argv[3]);
    printf("dimension K=%d  number of entries M=%d\n", k, m);
    vec = (COMP*)malloc(sizeof(COMP)*k);
    cb = (COMP*)malloc(sizeof(COMP)*k*m);
    cent = (COMP*)malloc(sizeof(COMP)*k*m);
    n = (int*)malloc(sizeof(int)*m);
    if (cb == NULL || cb == NULL || cent == NULL || vec == NULL) {
	printf("Error in malloc.\n");
	exit(1);
    }

    /* determine size of training set */

    J = 0;
    while(fread(vec, sizeof(COMP), k, ftrain) == (size_t)k)
	J++;
    printf("J=%d entries in training set\n", J);

    /* set up initial codebook state from samples of training set */

    rewind(ftrain);
    ret = fread(cb, sizeof(COMP), k*m, ftrain);
    //print_vec(cb, k, m);

    /* main loop */

    Dn = 1E32;
    j = 1;
    do {
	Dn_1 = Dn;

	/* zero centroids */

	for(i=0; i<m; i++) {
	    zero(&cent[i*k], k);
	    n[i] = 0;
	}

	/* quantise training set */

	se = 0.0;
	rewind(ftrain);
	for(i=0; i<J; i++) {
	    ret = fread(vec, sizeof(COMP), k, ftrain);
	    ind = quantise(cb, vec, k, m, &se);
	    //printf("ind %d\n", ind);
	    n[ind]++;
	    acc(&cent[ind*k], vec, k);
	}
	Dn = se/J;
	delta = (Dn_1-Dn)/Dn;

	printf("\r  Iteration %d, Dn = %f, Delta = %e\n", j, Dn, delta);
	j++;

	/* determine new codebook from centroids */

	for(i=0; i<m; i++) {
	    norm(&cent[i*k], k, n[i]);
	    memcpy(&cb[i*k], &cent[i*k], k*sizeof(COMP));
	}

    } while (delta > DELTAQ);

    /* save codebook to disk */

    fvq = fopen(argv[4],"wt");
    if (fvq == NULL) {
	printf("Error opening VQ file: %s\n",argv[4]);
	exit(1);
    }

    fprintf(fvq,"%d %d\n",k,m);
    for(j=0; j<m; j++) {
	for(i=0; i<k; i++)
	    fprintf(fvq,"%f %f ", cb[j*k+i].real, cb[j*k+i].imag);
	fprintf(fvq,"\n");
    }
    fclose(fvq);

    return 0;
}

/*-----------------------------------------------------------------------*\

				FUNCTIONS

\*-----------------------------------------------------------------------*/

void print_vec(COMP cb[], int k, int m)
{
    int i,j;

    for(j=0; j<m; j++) {
	for(i=0; i<k; i++) 
	    printf("%f %f ", cb[j*k+i].real, cb[j*k+i].imag);
	printf("\n");
    }
}


static COMP cconj(COMP a)
{
    COMP res;

    res.real = a.real;
    res.imag = -a.imag;

    return res;
}

static COMP cmult(COMP a, COMP b)
{
    COMP res;

    res.real = a.real*b.real - a.imag*b.imag;
    res.imag = a.real*b.imag + a.imag*b.real;

    return res;
}

static COMP cadd(COMP a, COMP b)
{
    COMP res;

    res.real = a.real + b.real;
    res.imag = a.imag + b.imag;

    return res;
}

/*---------------------------------------------------------------------------*\

	FUNCTION....: zero()

	AUTHOR......: David Rowe
	DATE CREATED: 23/2/95

	Zeros a vector of length k.

\*---------------------------------------------------------------------------*/

void zero(COMP v[], int k)
{
    int	i;

    for(i=0; i<k; i++) {
	v[i].real = 0.0;
	v[i].imag = 0.0;
    }
}

/*---------------------------------------------------------------------------*\

	FUNCTION....: acc()

	AUTHOR......: David Rowe
	DATE CREATED: 23/2/95

	Adds k dimensional vectors v1 to v2 and stores the result back
	in v1.  We add them like vectors on the complex plane, summing
	the real and imag terms.

\*---------------------------------------------------------------------------*/

void acc(COMP v1[], COMP v2[], int k)
{
    int	   i;

    for(i=0; i<k; i++)
	v1[i] = cadd(v1[i], v2[i]);
}

/*---------------------------------------------------------------------------*\

	FUNCTION....: norm()

	AUTHOR......: David Rowe
	DATE CREATED: 23/2/95

	Normalises each element in k dimensional vector.

\*---------------------------------------------------------------------------*/

void norm(COMP v[], int k, int n)
{
    int	   i;
    float  mag;

    for(i=0; i<k; i++) {
	mag = sqrt(v[i].real*v[i].real + v[i].imag*v[i].imag);
	v[i].real /= mag;
	v[i].imag /= mag;
    }
}

/*---------------------------------------------------------------------------*\

	FUNCTION....: quantise()

	AUTHOR......: David Rowe
	DATE CREATED: 23/2/95

	Quantises vec by choosing the nearest vector in codebook cb, and
	returns the vector index.  The squared error of the quantised vector
	is added to se.

\*---------------------------------------------------------------------------*/

int quantise(COMP cb[], COMP vec[], int k, int m, float *se)
{
   float   e;		/* current error		*/
   long	   besti;	/* best index so far		*/
   float   beste;	/* best error so far		*/
   long	   j;
   int     i;
   COMP    diff;

   besti = 0;
   beste = 1E32;
   for(j=0; j<m; j++) {
	e = 0.0;
	for(i=0; i<k; i++) {
	    diff = cmult(cb[j*k+i], cconj(vec[i]));
	    e += pow(atan2(diff.imag, diff.real), 2.0);
	}
	if (e < beste) {
	    beste = e;
	    besti = j;
	}
   }

   *se += beste;

   return(besti);
}

