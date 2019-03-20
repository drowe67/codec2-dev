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

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <getopt.h>

/*-----------------------------------------------------------------------*\

				DEFINES

\*-----------------------------------------------------------------------*/

#define	DELTAQ 	0.005		/* quiting distortion			*/
#define	MAX_STR	80		/* maximum string length		*/

/*-----------------------------------------------------------------------*\

			FUNCTION PROTOTYPES

\*-----------------------------------------------------------------------*/

void zero(float v[], int k);
void acc(float v1[], float v2[], int k);
void norm(float v[], int k, long n);
long quantise(float cb[], float vec[], int k, int m, float *beste, float *se);

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
    float  e;           /* sqaured error for current vector             */
    float  se;		/* squared error for this iteration		*/
    float  var,var_1;	/* current and previous iterations distortion	*/
    float  delta;	/* improvement in distortion 			*/
    long   noutliers[3];/* number of vectors quantisers with > 3*sd     */
    FILE   *ftrain;	/* file containing training set			*/
    FILE   *fvq;	/* file containing vector quantiser		*/
    int     ret;
    float   deltaq_stop = DELTAQ;
    FILE   *fres = NULL;
    
    int o = 0;
    int opt_idx = 0;
    while( o != -1 ) {
        static struct option long_opts[] = {
            {"help",     no_argument,       0, 'h'},
            {"residual", required_argument, 0, 'r'},
            {"stop",     required_argument, 0, 's'},
            {0, 0, 0, 0}
        };
        
        o = getopt_long(argc,argv,"hr:s:",long_opts,&opt_idx);
        
        switch(o) {
        case 'r':
            fres = fopen(optarg,"wb"); assert(fres != NULL);
            //fprintf(stderr, "writing res to : %s \n", optarg);
            break;
        case 's':
            deltaq_stop = atof(optarg);
            //fprintf(stderr, "deltaq_stop :%f\n", deltaq_stop);
            break;
        case 'h':
        case '?':
            goto helpmsg;
            break;
        }
    }
    int dx = optind;

    //fprintf(stderr, "argc: %d dx: %d\n", argc, dx);
    if ((argc - dx) < 4) {
        fprintf(stderr, "Too few arguments\n");
    helpmsg:
        fprintf(stderr, "usage: %s [Options] TrainFile.f32 K(dimension) M(codebook size) VQFile.f32\n", argv[0]);
        fprintf(stderr, "  -r --residual VQResidualErrorFile.f32usage\n");
        fprintf(stderr, "  -s --stop StopDelta\n");
        exit(1);
    }

    /* Open training file */

    ftrain = fopen(argv[dx],"rb");
    if (ftrain == NULL) {
	printf("Error opening training database file: %s\n",argv[dx]);
	exit(1);
    }

    /* determine k and m, and allocate arrays */

    k = atol(argv[dx+1]);
    m = atol(argv[dx+2]);
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
        quantise(cb, vec, k, 1, &e, &se);
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

	se = 0.0; noutliers[0] = noutliers[1] = noutliers[2] = 0;
	rewind(ftrain);
	for(i=0; i<J; i++) {
	    ret = fread(vec, sizeof(float), k, ftrain);
            assert(ret == k);
	    ind = quantise(cb, vec, k, m, &e, &se);
	    n[ind]++;
	    acc(&cent[ind*k], vec, k);
            //if (i < 100)
            //    printf("e: %f sqrt(e/k): %f sd: %f noutliers: %ld\n", e, sqrt(e/k), sd, noutliers[0]);
            if (sqrt(e/k) > 1.0) noutliers[0]++;
            if (sqrt(e/k) > 2.0) noutliers[1]++;
            if (sqrt(e/k) > 3.0) noutliers[2]++;
	}
	var = se/(J*k);
	delta = (var_1-var)/var;

	printf("\r  Iteration %ld, var = %4.2f, sd = %4.2f outliers > 1/2/3 dB = %3.2f/%f3.2/%3.2f Delta = %5.4f\n", j, var, sqrt(var),
               (float)noutliers[0]/J, (float)noutliers[1]/J, (float)noutliers[2]/J, delta);
	j++;

	/* determine new codebook from centroids */

	if (delta > deltaq_stop)
	    for(i=0; i<m; i++) {
		if (n[i] != 0) {
		    norm(&cent[i*k], k, n[i]);
		    memcpy(&cb[i*k], &cent[i*k], k*sizeof(float));
		}
	    }

    } while (delta > deltaq_stop);

    /* save VQ to disk */

    fvq = fopen(argv[dx+3],"wt");
    if (fvq == NULL) {
	printf("Error opening VQ file: %s\n",argv[dx+3]);
	exit(1);
    }

    fwrite(cb, sizeof(float), m*k, fvq);
    
    /* optionally output residual error for next stage VQ */

    if (fres != NULL) {
        float res[k];
	rewind(ftrain);
	for(j=0; j<J; j++) {
	    ret = fread(vec, sizeof(float), k, ftrain);
	    ind = quantise(cb, vec, k, m, &e, &se);
            for(i=0; i<k; i++) {
                res[i] = vec[i] - cb[k*ind+i];
            }
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

long quantise(float cb[], float vec[], int k, int m, float *beste, float *se)
/* float   cb[][K];	current VQ codebook		*/
/* float   vec[];	vector to quantise		*/
/* int	   k;		dimension of vectors		*/
/* int     m;		size of codebook		*/
/* float   beste;	current squared error		*/
/* float   *se;		accumulated squared error 	*/
{
   long	   besti;	/* best index so far		*/
   long	   j;
   int     i;
   float   diff,e;

   besti = 0;
   *beste = 1E32;
   for(j=0; j<m; j++) {
	e = 0.0;
	for(i=0; i<k; i++) {
	    diff = cb[j*k+i]-vec[i];
	    e += diff*diff;
	}
	if (e < *beste) {
	    *beste = e;
	    besti = j;
	}
   }

   *se += *beste;
   
   return(besti);
}

