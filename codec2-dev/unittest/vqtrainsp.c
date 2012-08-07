/*--------------------------------------------------------------------------*\

	FILE........: vqtrainsp.c
	AUTHOR......: David Rowe
	DATE CREATED: 7 August 2012

	This program trains sparse amplitude vector quantisers.
	Modified from vqtrainph.c

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

void zero(float v[], int d);
void acc(float v1[], float v2[], int d);
void norm(float v[], int k, int n[]);
int quantise(float cb[], float vec[], int d, int e, float *se);
void print_vec(float cb[], int d, int e);

/*-----------------------------------------------------------------------* \

				MAIN

\*-----------------------------------------------------------------------*/

int main(int argc, char *argv[]) {
    int    d,e;		/* dimension and codebook size			*/
    float  *vec;	/* current vector 				*/
    float  *cb;		/* vector codebook				*/
    float  *cent;	/* centroids for each codebook entry		*/
    int    *n;		/* number of vectors in this interval		*/
    int     J;		/* number of vectors in training set		*/
    int     ind;	/* index of current vector			*/
    float   se;	        /* total squared error for this iteration       */
    float   var;        /* variance                                     */ 
    float   var_1;	/* previous variance            	        */
    float   delta;	/* improvement in distortion 			*/
    FILE   *ftrain;	/* file containing training set			*/
    FILE   *fvq;	/* file containing vector quantiser		*/
    int     ret;
    int     i,j, finished, iterations;
    float   sd;
    int     var_n;

    /* Interpret command line arguments */

    if (argc != 5)	{
	printf("usage: %s TrainFile D(dimension) E(number of entries) VQFile\n", argv[0]);
	exit(1);
    }

    /* Open training file */

    ftrain = fopen(argv[1],"rb");
    if (ftrain == NULL) {
	printf("Error opening training database file: %s\n",argv[1]);
	exit(1);
    }

    /* determine k and m, and allocate arrays */

    d = atoi(argv[2]);
    e = atoi(argv[3]);
    printf("\n");
    printf("dimension D=%d  number of entries E=%d\n", d, e);
    vec = (float*)malloc(sizeof(float)*d);
    cb = (float*)malloc(sizeof(float)*d*e);
    cent = (float*)malloc(sizeof(float)*d*e);
    n = (int*)malloc(sizeof(int)*d*e);
    if (cb == NULL || cb == NULL || cent == NULL || vec == NULL) {
	printf("Error in malloc.\n");
	exit(1);
    }

    /* determine size of training set */

    J = 0;
    var_n = 0;
    while(fread(vec, sizeof(float), d, ftrain) == (size_t)d) {
	for(j=0; j<d; j++)
	    if (vec[j] != 0.0)
		var_n++;
	J++;
    }
    printf("J=%d sparse vectors in training set, %d non-zero values\n", J, var_n);

    /* set up initial codebook state from samples of training set */

    #define INIT1
    #ifdef INIT1
    rewind(ftrain);
    ret = fread(cb, sizeof(float), d*e, ftrain);
    #endif

    #ifdef INIT2
    for(i=0; i<d*e; i++)
	cb[i] = 1.0 - 2.0*rand()/RAND_MAX;
    #endif

    //print_vec(cb, d, 2);

    /* main loop */

    printf("\n");
    printf("Iteration  delta  var    std dev\n");
    printf("--------------------------------\n");

    iterations = 0;
    finished = 0;
    delta = 0;
    var_1 = 0.0;

    do {
	/* zero centroids */

	for(i=0; i<e; i++) {
	    zero(&cent[i*d], d);
	    for(j=0; j<d; j++)
		n[i*d+j] = 0;
	}

	//#define DBG
	#ifdef DBG
	printf("cb...\n");
	print_vec(cb, d, e);
	printf("\n\nquantise...\n");
	#endif

	/* quantise training set */

	se = 0.0;
	rewind(ftrain);
	for(i=0; i<J; i++) {
	    ret = fread(vec, sizeof(float), d, ftrain);
	    ind = quantise(cb, vec, d, e, &se);
	    #ifdef DBG
	    print_vec(vec, d, 1);
	    printf("      ind %d se: %f\n", ind, se);
	    #endif
	    acc(&cent[ind*d], vec, d);
	    for(j=0; j<d; j++)
		if (vec[j] != 0.0)
		    n[ind*d+j]++;
	}
	
	#ifdef DBG
	printf("cent...\n");
	print_vec(cent, d, e);
	printf("\n");
	#endif

	/* work out stats */

	var = se/var_n;	
	sd = sqrt(var);

	iterations++;
	if (iterations > 1) {
	    if (var > 0.0) {
		delta = (var_1 - var)/var;
	    }
	    else
		delta = 0;
	    if (delta < DELTAQ)
		finished = 1;
	}      
		     
	/* determine new codebook from centroids */

	for(i=0; i<e; i++) {
	    norm(&cent[i*d], d, &n[i*d]);
	    memcpy(&cb[i*d], &cent[i*d], d*sizeof(float));
	}
	
	#ifdef DBG
	printf("new cb ...\n");
	print_vec(cent, d, e);
	printf("\n");
	#endif

	printf("%2d         %4.3f  %4.3f  %4.3f \n",iterations, delta, var, sd);
	
	var_1 = var;
    } while (!finished);


    //print_vec(cb, d, 1);
    
    /* save codebook to disk */

    fvq = fopen(argv[4],"wt");
    if (fvq == NULL) {
	printf("Error opening VQ file: %s\n",argv[4]);
	exit(1);
    }

    fprintf(fvq,"%d %d\n",d,e);
    for(j=0; j<e; j++) {
	for(i=0; i<d; i++)
	    fprintf(fvq,"% 4.3f ", cb[j*d+i]);
	fprintf(fvq,"\n");
    }
    fclose(fvq);

    return 0;
}

/*-----------------------------------------------------------------------*\

				FUNCTIONS

\*-----------------------------------------------------------------------*/

void print_vec(float cb[], int d, int e)
{
    int i,j;

    for(j=0; j<e; j++) {
	printf("    ");
	for(i=0; i<d; i++) 
	    printf("% 4.3f ", cb[j*d+i]);
	printf("\n");
    }
}


/*---------------------------------------------------------------------------*\

	FUNCTION....: zero()

	AUTHOR......: David Rowe
	DATE CREATED: 23/2/95

	Zeros a vector of length d.

\*---------------------------------------------------------------------------*/

void zero(float v[], int d)
{
    int	i;

    for(i=0; i<d; i++) {
	v[i] = 0.0;
    }
}

/*---------------------------------------------------------------------------*\

	FUNCTION....: acc()

	AUTHOR......: David Rowe
	DATE CREATED: 23/2/95

	Adds d dimensional vectors v1 to v2 and stores the result back
	in v1.  

	An unused entry in a sparse vector is set to zero so won't
	affect the accumulation process.

\*---------------------------------------------------------------------------*/

void acc(float v1[], float v2[], int d)
{
    int	   i;

    for(i=0; i<d; i++)
	v1[i] += v2[i];
}

/*---------------------------------------------------------------------------*\

	FUNCTION....: norm()

	AUTHOR......: David Rowe
	DATE CREATED: 23/2/95

	Normalises each element in d dimensional vector.

\*---------------------------------------------------------------------------*/

void norm(float v[], int d, int n[])
{
    int	   i;

    for(i=0; i<d; i++) {
	if (n[i] != 0)
	    v[i] /= n[i];
    }
}

/*---------------------------------------------------------------------------*\

	FUNCTION....: quantise()

	AUTHOR......: David Rowe
	DATE CREATED: 23/2/95

	Quantises vec by choosing the nearest vector in codebook cb, and
	returns the vector index.  The squared error of the quantised vector
	is added to se.  

	Unused entries in sparse vectors are ignored.

\*---------------------------------------------------------------------------*/

int quantise(float cb[], float vec[], int d, int e, float *se)
{
   float   error;	/* current error		*/
   int     besti;	/* best index so far		*/
   float   best_error;	/* best error so far		*/
   int	   i,j,n;
   float   diff;

   besti = 0;
   best_error = 1E32;
   for(j=0; j<e; j++) {
       error = 0.0; n = 0;
       for(i=0; i<d; i++) {
	   if ((vec[i] != 0.0) && (cb[j*d+i] != 0.0)) {
	       diff = cb[j*d+i] - vec[i];
	       error += diff*diff;
	       n++;
	   }
       }
       error /= n;
       if (error < best_error) {
	   best_error = error;
	   besti = j;
       }
   }

   *se += best_error;

   return(besti);
}

