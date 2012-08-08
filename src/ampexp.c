/*---------------------------------------------------------------------------*\
                                                                             
  FILE........: ampexp.c                                
  AUTHOR......: David Rowe                                             
  DATE CREATED: 7 August 2012
                                                                             
  Functions for experimenting with amplitude quantisation.
                                                                             
\*---------------------------------------------------------------------------*/

/*
  Copyright (C) 2012 David Rowe

  All rights reserved.

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU Lesser General Public License version 2.1, as
  published by the Free Software Foundation.  This program is
  distributed in the hope that it will be useful, but WITHOUT ANY
  WARRANTY; without even the implied warranty of MERCHANTABILITY or
  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public
  License for more details.

  You should have received a copy of the GNU Lesser General Public License
  along with this program; if not,see <http://www.gnu.org/licenses/>. 
*/


#include <assert.h>
#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "ampexp.h"


#define PRED_COEFF 0.9

/* states for amplitude experiments */

struct codebook {
    unsigned int	 k;
    unsigned int	 log2m;
    unsigned int	 m;
    float               *cb;
    unsigned int         offset; 
};

struct AEXP {
    float            A_prev[MAX_AMP];
    int              frames;
    float            snr;
    int              snr_n;
    float            var;
    int              var_n;
    float            vq_var;
    int              vq_var_n;
    struct codebook *vq;
    int              offset;
};


/*---------------------------------------------------------------------------*\

  Bruce Perens' funcs to load codebook files

\*---------------------------------------------------------------------------*/


static const char format[] =
"The table format must be:\n"
"\tTwo integers describing the dimensions of the codebook.\n"
"\tThen, enough numbers to fill the specified dimensions.\n";

static float get_float(FILE * in, const char * name, char * * cursor, char * buffer, int size)
{
  for ( ; ; ) {
    char *	s = *cursor;
    char	c;

    while ( (c = *s) != '\0' && !isdigit(c) && c != '-' && c != '.' )
      s++;
     
    /* Comments start with "#" and continue to the end of the line. */
    if ( c != '\0' && c != '#' ) {
      char *	end = 0;
      float	f = 0;

      f = strtod(s, &end);

      if ( end != s )
        *cursor = end;
        return f;
    }

    if ( fgets(buffer, size, in) == NULL ) {
      fprintf(stderr, "%s: Format error. %s\n", name, format);
      exit(1);
    }
    *cursor = buffer;
  }
}

static struct codebook *load(const char * name)
{
    FILE               *file;
    char		line[2048];
    char               *cursor = line;
    struct codebook    *b = malloc(sizeof(struct codebook));
    int			i;
    int			size;

    file = fopen(name, "rt");
    assert(file != NULL);

    *cursor = '\0';

    b->k = (int)get_float(file, name, &cursor, line, sizeof(line));
    b->m = (int)get_float(file, name ,&cursor, line, sizeof(line));
    size = b->k * b->m;

    b->cb = (float *)malloc(size * sizeof(float));

    for ( i = 0; i < size; i++ ) {
	b->cb[i] = get_float(file, name, &cursor, line, sizeof(line));
    }

    fclose(file);

    return b;
}


/*---------------------------------------------------------------------------* \

  amp_experiment_create()

  Inits states for amplitude quantisation experiments.

\*---------------------------------------------------------------------------*/

struct AEXP *amp_experiment_create() {
    struct AEXP *aexp;
    int i;

    aexp = (struct AEXP *)malloc(sizeof(struct AEXP));
    assert (aexp != NULL);

    for(i=0; i<MAX_AMP; i++)
	aexp->A_prev[i] = 1.0;
    aexp->frames = 0;
    aexp->snr = 0.0;
    aexp->snr_n = 0;
    aexp->var = 0.0;
    aexp->var_n = 0;
    aexp->vq_var = 0.0;
    aexp->vq_var_n = 0;
    aexp->vq = load("../unittest/amp1_20_1024.txt");
    //aexp->vq = load("hts1a_vq1_20_300.txt");
    aexp->offset = 20;

    return aexp;
}


/*---------------------------------------------------------------------------* \

  amp_experiment_destroy()

\*---------------------------------------------------------------------------*/

void amp_experiment_destroy(struct AEXP *aexp) {
    assert(aexp != NULL);
    if (aexp->snr != 0.0)
	printf("snr: %4.2f dB\n", aexp->snr/aexp->snr_n);
    if (aexp->var != 0.0)
	printf("var...: %4.3f  std dev...: %4.3f (%d amplitude samples)\n", 
	       aexp->var/aexp->var_n, sqrt(aexp->var/aexp->var_n), aexp->var_n);
    if (aexp->vq_var != 0.0)
	printf("vq var: %4.3f  std dev...: %4.3f (%d amplitude samples)\n", 
	       aexp->vq_var/aexp->vq_var_n, sqrt(aexp->vq_var/aexp->vq_var_n), aexp->vq_var_n);
    free(aexp);
}


/*---------------------------------------------------------------------------*\

  Various test and experimental functions ................

\*---------------------------------------------------------------------------*/

/*
  Quantisation noise simulation.  Assume noise on amplitudes is a uniform
  distribution, of +/- x dB.  This means x = sqrt(3)*sigma.

  Note: for uniform distribution var = = sigma * sigma = (b-a)*(b-a)/12. 
*/

static void add_quant_noise(struct AEXP *aexp, MODEL *model, int start, int end, float sigma_dB)
{
    int   m;
    float x_dB;
    float noise_sam_dB;
    float noise_sam_lin;

    x_dB = sqrt(3.0) * sigma_dB;

    for(m=start; m<=end; m++) {
	noise_sam_dB = x_dB*(1.0 - 2.0*rand()/RAND_MAX);
	//printf("%f\n", noise_sam_dB);
	noise_sam_lin = pow(10.0, noise_sam_dB/20.0);
	model->A[m] *= noise_sam_lin;
	aexp->var += noise_sam_dB*noise_sam_dB;
	aexp->var_n++;
    }

}

/*
  void print_sparse_pred_error()

  use to check pred error stats (e.g. of first 1kHz) in Octave:

     $ ./c2sim ../raw/hts1a.raw --ampexp > amppe.txt

     octave> load ../src/amppe.txt
     octave> std(nonzeros(amppe(:,1:20)))
     octave> hist(nonzeros(amppe(:,1:20)),20);

 */


static void print_sparse_pred_error(struct AEXP *aexp, MODEL *model, float mag_thresh)
{
    int    m, index;
    float  mag, error;
    float  sparse_pe[MAX_AMP];

    mag = 0.0;
    for(m=1; m<=model->L; m++)
	mag += model->A[m]*model->A[m];
    mag = 10*log10(mag/model->L);
    
    if (mag > mag_thresh) {
	for(m=0; m<MAX_AMP; m++) {
	    sparse_pe[m] = 0.0;
	}

	for(m=1; m<=model->L; m++) {
	    assert(model->A[m] > 0.0);
	    error = PRED_COEFF*20.0*log10(aexp->A_prev[m]) - 20.0*log10(model->A[m]);

	    index = MAX_AMP*m*model->Wo/PI;
	    assert(index < MAX_AMP);
	    sparse_pe[index] = error;
	}

	/* dump spare phase vector */

	for(m=0; m<MAX_AMP/4; m++)
	    printf("%f ", sparse_pe[m]);
	printf("\n");
    }
}


int vq_amp(float cb[], float vec[], float weights[], int d, int e, float *se)
{
   float   error;	/* current error		*/
   int     besti;	/* best index so far		*/
   float   best_error;	/* best error so far		*/
   int	   i,j;
   float   diff, metric, best_metric;

   besti = 0;
   best_metric = 1E32;
   for(j=0; j<e; j++) {
       metric = error = 0.0;
       for(i=0; i<d; i++) {
	   if (vec[i] != 0.0) {
	       diff = (cb[j*d+i] - vec[i]);
	       error += diff*diff;
	       metric += weights[i]*diff*diff;
	   }
       }
       if (metric < best_metric) {
	   best_error = error;
	   best_metric = metric;
	   besti = j;
       }
   }

   *se += best_error;

   return(besti);
}


static void split_vq(float sparse_pe_out[], struct AEXP *pexp, struct codebook *vq, float weights[], float sparse_pe_in[])
{
    int i, j, non_zero, vq_ind;
    float se;

    vq_ind = vq_amp(vq->cb, &sparse_pe_in[vq->offset], &weights[vq->offset], vq->k, vq->m, &se);
    printf("\n offset %d k %d m %d vq_ind %d j: ", vq->offset, vq->k, vq->m, vq_ind);
  
    non_zero = 0;
    for(i=0, j=vq->offset; i<vq->k; i++,j++) {
	if (sparse_pe_in[j] != 0.0) {
	    //printf("%d ", j);
	    sparse_pe_out[j] = vq->cb[vq->k * vq_ind + i];
	    non_zero++;
	}
    }
    pexp->vq_var_n += non_zero;
}


static void sparse_vq_pred_error(struct AEXP *aexp, 
				 MODEL       *model 
)
{
    int    m, index;
    float  error, amp_dB, mag;
    float  sparse_pe_in[MAX_AMP];
    float  sparse_pe_out[MAX_AMP];
    float  weights[MAX_AMP];

    for(m=0; m<MAX_AMP; m++) {
	sparse_pe_in[m] = 0.0;
	sparse_pe_out[m] = 0.0;
    }

    mag = 0.0;
    for(m=1; m<=model->L; m++) {
	mag += model->A[m]*model->A[m];
 	assert(model->A[m] > 0.0);
	error = PRED_COEFF*20.0*log10(aexp->A_prev[m]) - 20.0*log10(model->A[m]);

	index = MAX_AMP*m*model->Wo/PI;
	assert(index < MAX_AMP);
	sparse_pe_in[index] = error;
	weights[index] = model->A[m];
    }
    mag = 10.0*log10(mag/model->L);

    /* vector quantise */
        
    for(m=0; m<MAX_AMP; m++) {
	sparse_pe_out[m] = sparse_pe_in[m];
    }

    split_vq(sparse_pe_out, aexp, aexp->vq, weights, sparse_pe_in);
    #ifdef SIM_VQ
    for(m=0; m<MAX_AMP; m++) {
	if (sparse_pe_in[m] != 0.0) {
	    float error = 4.0*(1.0 - 2.0*rand()/RAND_MAX);
	    aexp->vq_var += error*error;
	    aexp->vq_var_n++;
	    sparse_pe_out[m] = sparse_pe_in[m] + error;
	}
    }
    #endif

    if (mag > 40.0)
	for(m=0; m<MAX_AMP; m++) {
	    if (sparse_pe_in[m] != 0.0)
		aexp->vq_var += pow(sparse_pe_out[m] - sparse_pe_in[m], 2.0);
	}
    
    /* transform quantised amps back */

    for(m=1; m<=model->L; m++) {
	index = MAX_AMP*m*model->Wo/PI;
	assert(index < MAX_AMP);
	amp_dB = PRED_COEFF*20.0*log10(aexp->A_prev[m]) - sparse_pe_out[index];
	//printf("in: %f  out: %f\n", sparse_pe_in[index], sparse_pe_out[index]);
	//printf("amp_dB: %f A[m] (dB) %f\n", amp_dB, 20.0*log10(model->A[m]));
	model->A[m] = pow(10.0, amp_dB/20.0);
    }
    //exit(0);
}


static void update_snr_calc(struct AEXP *aexp, MODEL *model, float before[])
{
    int m;
    float signal, noise, signal_dB;

    signal = 0.0; noise = 0.0;
    for(m=1; m<=model->L; m++) {	    
	signal += before[m]*before[m];
	noise  += pow(before[m] - model->A[m], 2.0);
	//printf("%f %f\n", before[m], model->phi[m]);
    }
    signal_dB = 10*log10(signal/model->L);
    if (signal_dB > - 100.0) {
	aexp->snr += 10.0*log10(signal/noise);
	aexp->snr_n++;
    }
}


/*---------------------------------------------------------------------------*\

  amp_experiment()

  Amplitude quantisation experiments.

\*---------------------------------------------------------------------------*/

void amp_experiment(struct AEXP *aexp, MODEL *model) {
    int m;
    float before[MAX_AMP];

    for(m=1; m<=model->L; m++)
	before[m] = model->A[m];

    //print_sparse_pred_error(aexp, model, 40.0);
    //add_quant_noise(aexp, model, model->L/2, model->L, 3);
    sparse_vq_pred_error(aexp, model);

    update_snr_calc(aexp, model, before);

    /* update states */

    for(m=1; m<=model->L; m++)
	aexp->A_prev[m] = model->A[m];	    
    aexp->frames++;
}

