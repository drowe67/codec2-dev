/*--------------------------------------------------------------------------*\

	FILE........: VQTRAIN.C
	AUTHOR......: David Rowe
	DATE CREATED: 23/2/95

	This program trains vector quantisers using K dimensional Lloyd-Max
	method.

\*--------------------------------------------------------------------------*/

/*
  Copyright (C) 2009 David Rowe

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

/*-----------------------------------------------------------------------*\

				INCLUDES

\*-----------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <assert.h>

/*-----------------------------------------------------------------------*\

				DEFINES

\*-----------------------------------------------------------------------*/

#define	DELTAQ 	0.01		/* quiting distortion			*/
#define	MAX_STR	80		/* maximum string length		*/

/*-----------------------------------------------------------------------*\

			FUNCTION PROTOTYPES

\*-----------------------------------------------------------------------*/

void zero(float v[], int k);
void acc(float v1[], float v2[], int k);
void norm(float v[], int k, long n);
long quantise(float cb[], float vec[], int k, int m, float *se);

/*-----------------------------------------------------------------------* \

				MAIN

\*-----------------------------------------------------------------------*/

int main(int argc, char *argv[]) {
    long   k,m;		/* dimension and codebook size			*/
    float  *vec;	/* current vector 				*/
    float  *cb;		/* vector codebook				*/
    float  *cent;	/* centroids for each codebook entry		*/
    long   *n;		/* number of vectors in this interval		*/
    long   J;		/* number of vectors in training set		*/
    long   i,j;
    long   ind;	     	/* index of current vector			*/
    float  se;		/* squared error for this iteration		*/
    float  var,var_1;	/* current and previous iterations distortion	*/
    float  delta;	/* improvement in distortion 			*/
    FILE   *ftrain;	/* file containing training set			*/
    FILE   *fvq;	/* file containing vector quantiser		*/
    int     ret;

    /* Interpret command line arguments */

    if (argc < 5)	{
	printf("usage: %s TrainFile K(dimension) M(codebook size) VQFile [residual.f32]\n", argv[0]);
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
    printf("vector dimension K=%ld  codebook size M=%ld ", k, m);
    vec = (float*)malloc(sizeof(float)*k);
    cb = (float*)malloc(sizeof(float)*k*m);
    cent = (float*)malloc(sizeof(float)*k*m);
    n = (long*)malloc(sizeof(long)*m);
    if (cb == NULL || cb == NULL || cent == NULL || vec == NULL) {
	printf("Error in malloc.\n");
	exit(1);
    }

    /* determine size of training set */

    J = 0; zero(cent, k);
    while(fread(vec, sizeof(float), k, ftrain) == (size_t)k) {
        J++;
        acc(cent, vec, k);
    }
    printf("J=%ld vectors in training set\n", J);

    /* Interation is a 0 bit VQ (i.e. mean of training set) as starting point */
    
    norm(cent, k, J);
    memcpy(cb, cent, k*sizeof(float));
    se = 0.0;
    rewind(ftrain);
    for(i=0; i<J; i++) {
        ret = fread(vec, sizeof(float), k, ftrain);
        assert(ret == k);
        quantise(cb, vec, k, 1, &se);
    }
    var = se/(J*k);
    printf("\r  Iteration 0, var = %f, sd = %f\n", var, sqrt(var));

    /* set up initial codebook state from samples of training set */

    for(i=0; i<m; i++) {
        j = i*(J/m);
        fseek(ftrain, j*k*sizeof(float), SEEK_SET);
        ret = fread(&cb[i*k], sizeof(float), k, ftrain);
        assert(ret == k);
    }

    /* main loop */

    j = 1;
    do {
	var_1 = var;

	/* zero centroids */

	for(i=0; i<m; i++) {
	    zero(&cent[i*k], k);
	    n[i] = 0;
	}

	/* quantise training set */

	se = 0.0;
	rewind(ftrain);
	for(i=0; i<J; i++) {
	    ret = fread(vec, sizeof(float), k, ftrain);
            assert(ret == k);
	    ind = quantise(cb, vec, k, m, &se);
	    n[ind]++;
	    acc(&cent[ind*k], vec, k);
	}
	var = se/(J*k);
	delta = (var_1-var)/var;

	printf("\r  Iteration %ld, var = %f, sd = %f Delta = %e\n", j, var, sqrt(var), delta);
	j++;

	/* determine new codebook from centroids */

	if (delta > DELTAQ)
	    for(i=0; i<m; i++) {
		if (n[i] != 0) {
		    norm(&cent[i*k], k, n[i]);
		    memcpy(&cb[i*k], &cent[i*k], k*sizeof(float));
		}
	    }

    } while (delta > DELTAQ);

    /* save codebook to disk */

    fvq = fopen(argv[4],"wt");
    if (fvq == NULL) {
	printf("Error opening VQ file: %s\n",argv[4]);
	exit(1);
    }

    fprintf(fvq,"%ld %ld\n",k,m);
    for(j=0; j<m; j++) {
	for(i=0; i<k; i++)
	    fprintf(fvq,"%f  ",cb[j*k+i]);
	fprintf(fvq,"\n");
    }

    /* optionally output for next stage VQ */
    
    if (argc == 6) {
        FILE *fres = fopen(argv[5],"wb"); assert(fres != NULL);
        float res[k];
	rewind(ftrain);
	for(j=0; j<J; j++) {
	    ret = fread(vec, sizeof(float), k, ftrain);
	    ind = quantise(cb, vec, k, m, &se);
            for(i=0; i<k; i++) {
                res[i] = vec[i] - cb[k*ind+i];
                if (j<2) printf("%3.2f ", res[i]);
            }
            if (j<2) printf("\n");
            fwrite(res, sizeof(float), k, fres);
	}
        fclose(fres);
    }
    
    fclose(fvq);
    fclose(ftrain);
    free(vec);
    free(n);

    return 0;
}

/*-----------------------------------------------------------------------*\

				FUNCTIONS

\*-----------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*\

	FUNCTION....: zero()

	AUTHOR......: David Rowe
	DATE CREATED: 23/2/95

	Zeros a vector of length k.

\*---------------------------------------------------------------------------*/

void zero(float v[], int k)
/*  float  v[];		ptr to start of vector		*/
/*  int    k;		lngth of vector			*/
{
    int	i;

    for(i=0; i<k; i++)
	v[i] = 0.0;
}

/*---------------------------------------------------------------------------*\

	FUNCTION....: acc()

	AUTHOR......: David Rowe
	DATE CREATED: 23/2/95

	Adds k dimensional vectors v1 to v2 and stores the result back in v1.

\*---------------------------------------------------------------------------*/

void acc(float v1[], float v2[], int k)
/*  float  v1[];	ptr to start of vector to accumulate	*/
/*  float  v2[];	ptr to start of vector to add		*/
/*  int	   k;		dimension of vectors			*/
{
    int	   i;

    for(i=0; i<k; i++)
	v1[i] += v2[i];
}

/*---------------------------------------------------------------------------*\

	FUNCTION....: norm()

	AUTHOR......: David Rowe
	DATE CREATED: 23/2/95

	Divides each element in k dimensional vector v by n.

\*---------------------------------------------------------------------------*/

void norm(float v[], int k, long n)
/*  float  v[];		ptr to start of vector		*/
/*  int	   k;		dimension of vectors		*/
/*  int	   n;		normalising factor		*/
{
    int	   i;

    assert(n != 0);
    for(i=0; i<k; i++)
	v[i] /= n;
}

/*---------------------------------------------------------------------------*\

	FUNCTION....: quantise()

	AUTHOR......: David Rowe
	DATE CREATED: 23/2/95

	Quantises vec by choosing the nearest vector in codebook cb, and
	returns the vector index.  The squared error of the quantised vector
	is added to se.

\*---------------------------------------------------------------------------*/

long quantise(float cb[], float vec[], int k, int m, float *se)
/* float   cb[][K];	current VQ codebook		*/
/* float   vec[];	vector to quantise		*/
/* int	   k;		dimension of vectors		*/
/* int     m;		size of codebook		*/
/* float   *se;		accumulated squared error 	*/
{
   float   e;		/* current error		*/
   long	   besti;	/* best index so far		*/
   float   beste;	/* best error so far		*/
   long	   j;
   int     i;
   float   diff;

   besti = 0;
   beste = 1E32;
   for(j=0; j<m; j++) {
	e = 0.0;
	for(i=0; i<k; i++) {
	    diff = cb[j*k+i]-vec[i];
	    e += diff*diff;
	}
	if (e < beste) {
	    beste = e;
	    besti = j;
	}
   }

   *se += beste;

   return(besti);
}

