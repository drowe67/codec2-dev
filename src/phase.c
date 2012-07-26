/*---------------------------------------------------------------------------*\
                                                                             
  FILE........: phase.c                                           
  AUTHOR......: David Rowe                                             
  DATE CREATED: 1/2/09                                                 
                                                                             
  Functions for modelling and synthesising phase.
                                                                             
\*---------------------------------------------------------------------------*/

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
  along with this program; if not,see <http://www.gnu.org/licenses/>. 
*/

#include "defines.h"
#include "phase.h"
#include "kiss_fft.h"
#include "comp.h"
#include "glottal.c"

#include <assert.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

#define GLOTTAL_FFT_SIZE 512

/*---------------------------------------------------------------------------*\

  aks_to_H()

  Samples the complex LPC synthesis filter spectrum at the harmonic
  frequencies.

\*---------------------------------------------------------------------------*/

void aks_to_H(
              kiss_fft_cfg fft_fwd_cfg, 
	      MODEL *model,	/* model parameters */
	      float  aks[],	/* LPC's */
	      float  G,	        /* energy term */
	      COMP   H[],	/* complex LPC spectral samples */
	      int    order
)
{
  COMP  pw[FFT_ENC];	/* power spectrum (input) */
  COMP  Pw[FFT_ENC];	/* power spectrum (output) */
  int   i,m;		/* loop variables */
  int   am,bm;		/* limits of current band */
  float r;		/* no. rads/bin */
  float Em;		/* energy in band */
  float Am;		/* spectral amplitude sample */
  int   b;		/* centre bin of harmonic */
  float phi_;		/* phase of LPC spectra */

  r = TWO_PI/(FFT_ENC);

  /* Determine DFT of A(exp(jw)) ------------------------------------------*/

  for(i=0; i<FFT_ENC; i++) {
    pw[i].real = 0.0;
    pw[i].imag = 0.0;
  }

  for(i=0; i<=order; i++)
    pw[i].real = aks[i];

  kiss_fft(fft_fwd_cfg, (kiss_fft_cpx *)pw, (kiss_fft_cpx *)Pw);

  /* Sample magnitude and phase at harmonics */

  for(m=1; m<=model->L; m++) {
    am = floor((m - 0.5)*model->Wo/r + 0.5);
    bm = floor((m + 0.5)*model->Wo/r + 0.5);
    b = floor(m*model->Wo/r + 0.5);

    Em = 0.0;
    for(i=am; i<bm; i++)
      Em += G/(Pw[i].real*Pw[i].real + Pw[i].imag*Pw[i].imag);
    Am = sqrt(fabs(Em/(bm-am)));

    phi_ = -atan2(Pw[b].imag,Pw[b].real);
    H[m].real = Am*cos(phi_);
    H[m].imag = Am*sin(phi_);
  }
}


/*---------------------------------------------------------------------------*\

   phase_synth_zero_order()

   Synthesises phases based on SNR and a rule based approach.  No phase 
   parameters are required apart from the SNR (which can be reduced to a
   1 bit V/UV decision per frame).

   The phase of each harmonic is modelled as the phase of a LPC
   synthesis filter excited by an impulse.  Unlike the first order
   model the position of the impulse is not transmitted, so we create
   an excitation pulse train using a rule based approach.  

   Consider a pulse train with a pulse starting time n=0, with pulses
   repeated at a rate of Wo, the fundamental frequency.  A pulse train
   in the time domain is equivalent to harmonics in the frequency
   domain.  We can make an excitation pulse train using a sum of
   sinsusoids:

     for(m=1; m<=L; m++)
       ex[n] = cos(m*Wo*n)

   Note: the Octave script ../octave/phase.m is an example of this if
   you would like to try making a pulse train.

   The phase of each excitation harmonic is:

     arg(E[m]) = mWo

   where E[m] are the complex excitation (freq domain) samples,
   arg(x), just returns the phase of a complex sample x.

   As we don't transmit the pulse position for this model, we need to
   synthesise it.  Now the excitation pulses occur at a rate of Wo.
   This means the phase of the first harmonic advances by N samples
   over a synthesis frame of N samples.  For example if Wo is pi/20
   (200 Hz), then over a 10ms frame (N=80 samples), the phase of the
   first harmonic would advance (pi/20)*80 = 4*pi or two complete
   cycles.

   We generate the excitation phase of the fundamental (first
   harmonic):

     arg[E[1]] = Wo*N;

   We then relate the phase of the m-th excitation harmonic to the
   phase of the fundamental as:

     arg(E[m]) = m*arg(E[1])

   This E[m] then gets passed through the LPC synthesis filter to
   determine the final harmonic phase.
     
   Comparing to speech synthesised using original phases:

   - Through headphones speech synthesised with this model is not as 
     good. Through a loudspeaker it is very close to original phases.

   - If there are voicing errors, the speech can sound clicky or
     staticy.  If V speech is mistakenly declared UV, this model tends to
     synthesise impulses or clicks, as there is usually very little shift or
     dispersion through the LPC filter.

   - When combined with LPC amplitude modelling there is an additional
     drop in quality.  I am not sure why, theory is interformant energy
     is raised making any phase errors more obvious.

   NOTES:

     1/ This synthesis model is effectively the same as a simple LPC-10
     vocoders, and yet sounds much better.  Why? Conventional wisdom
     (AMBE, MELP) says mixed voicing is required for high quality
     speech.

     2/ I am pretty sure the Lincoln Lab sinusoidal coding guys (like xMBE
     also from MIT) first described this zero phase model, I need to look
     up the paper.

     3/ Note that this approach could cause some discontinuities in
     the phase at the edge of synthesis frames, as no attempt is made
     to make sure that the phase tracks are continuous (the excitation
     phases are continuous, but not the final phases after filtering
     by the LPC spectra).  Technically this is a bad thing.  However
     this may actually be a good thing, disturbing the phase tracks a
     bit.  More research needed, e.g. test a synthesis model that adds
     a small delta-W to make phase tracks line up for voiced
     harmonics.

\*---------------------------------------------------------------------------*/

void phase_synth_zero_order(
    kiss_fft_cfg fft_fwd_cfg,     
    MODEL *model,
    float  aks[],
    float *ex_phase,            /* excitation phase of fundamental */
    int    order
)
{
  int   m;
  float new_phi;
  COMP  Ex[MAX_AMP];		/* excitation samples */
  COMP  A_[MAX_AMP];		/* synthesised harmonic samples */
  COMP  H[MAX_AMP];             /* LPC freq domain samples */
  float G;
  float jitter = 0.0;
  float r;
  int   b;

  G = 1.0;
  aks_to_H(fft_fwd_cfg, model, aks, G, H, order);

  /* 
     Update excitation fundamental phase track, this sets the position
     of each pitch pulse during voiced speech.  After much experiment
     I found that using just this frame's Wo improved quality for UV
     sounds compared to interpolating two frames Wo like this:
     
     ex_phase[0] += (*prev_Wo+mode->Wo)*N/2;
  */
  
  ex_phase[0] += (model->Wo)*N;
  ex_phase[0] -= TWO_PI*floor(ex_phase[0]/TWO_PI + 0.5);
  r = TWO_PI/GLOTTAL_FFT_SIZE;

  for(m=1; m<=model->L; m++) {

    /* generate excitation */

    if (model->voiced) {
	/* I think adding a little jitter helps improve low pitch
	   males like hts1a. This moves the onset of each harmonic
	   over at +/- 0.25 of a sample.
	*/
        jitter = 0.25*(1.0 - 2.0*rand()/RAND_MAX);
        b = floor(m*model->Wo/r + 0.5);
	if (b > ((GLOTTAL_FFT_SIZE/2)-1)) {
		b = (GLOTTAL_FFT_SIZE/2)-1;
	}
	Ex[m].real = cos(ex_phase[0]*m - jitter*model->Wo*m + glottal[b]);
	Ex[m].imag = sin(ex_phase[0]*m - jitter*model->Wo*m + glottal[b]);
    }
    else {

	/* When a few samples were tested I found that LPC filter
	   phase is not needed in the unvoiced case, but no harm in
	   keeping it.
        */
	float phi = TWO_PI*(float)rand()/RAND_MAX;
        Ex[m].real = cos(phi);
	Ex[m].imag = sin(phi);
    }

    /* filter using LPC filter */

    A_[m].real = H[m].real*Ex[m].real - H[m].imag*Ex[m].imag;
    A_[m].imag = H[m].imag*Ex[m].real + H[m].real*Ex[m].imag;

    /* modify sinusoidal phase */
   
    new_phi = atan2(A_[m].imag, A_[m].real+1E-12);
    model->phi[m] = new_phi;
  }

}

/* states for phase experiments */

struct PEXP {
    float phi_prev[MAX_AMP];
    float Wo_prev;
    int   frames;
};


/*---------------------------------------------------------------------------* \

  phase_experiment_create()

  Inits states for phase quantisation experiments.

\*---------------------------------------------------------------------------*/

struct PEXP * phase_experiment_create() {
    struct PEXP *pexp;
    int i;

    pexp = (struct PEXP *)malloc(sizeof(struct PEXP));
    assert (pexp != NULL);

    for(i=0; i<MAX_AMP; i++)
	pexp->phi_prev[i] = 0.0;
    pexp->Wo_prev = 0.0;
    pexp->frames = 0;

    return pexp;
}


/*---------------------------------------------------------------------------* \

  phase_experiment_destroy()

\*---------------------------------------------------------------------------*/

void phase_experiment_destroy(struct PEXP *pexp) {
    assert(pexp != NULL);
    free(pexp);
}


/*---------------------------------------------------------------------------* \

  Various test and experimental functions ................

\*---------------------------------------------------------------------------*/

/* Bubblesort to find highest amplitude harmonics */

struct AMPINDEX {
    float amp;
    int   index;
};

static void bubbleSort(struct AMPINDEX numbers[], int array_size)
{
    int i, j;
    struct AMPINDEX temp;
 
  for (i = (array_size - 1); i > 0; i--)
  {
    for (j = 1; j <= i; j++)
    {
      if (numbers[j-1].amp < numbers[j].amp)
      {
        temp = numbers[j-1];
        numbers[j-1] = numbers[j];
        numbers[j] = temp;
      }
    }
  }
}

    //#define AMP
    #ifdef AMP
    {
	int r = 0;
	float max = 0;
	for(i=1; i<model.L/4; i++)
	    if (model.A[i] > max) max = model.A[i];
	for(i=1; i<model.L/4; i++) {
	    if (model.A[i] < 0.25*max) {
		model.phi[i] += (PI/2)*(1.0 - 2.0*(float)rand()/RAND_MAX);
		r++;
	    }
	}
	printf("r %d L/4 %d\n", r, model.L/4);
    }
    #endif

   //#define AMP
    #ifdef AMP
    {
	int r = 0;
	float max = 0;
	for(i=1; i<model.L/4; i++)
	    if (model.A[i] > max) max = model.A[i];
	for(i=1; i<model.L/4; i++) {
	    if (model.A[i] < 0.25*max) {
		model.phi[i] += (PI/2)*(1.0 - 2.0*(float)rand()/RAND_MAX);
		r++;
	    }
	}
	printf("r %d L/4 %d\n", r, model.L/4);
    }
    #endif

    //#define AMP1
    #ifdef AMP1
    {
	int r = 0, nr = 0, st, dim;
	int top, j, max_harm, n_harm;
	struct AMPINDEX sorted[MAX_AMP];
	int ind[MAX_AMP];

	st = 1;
	dim = model.L/4;

	for(i=st,j=1; i<=st+dim; i++,j++) {
	    sorted[j].amp = model.A[i];
	    sorted[j].index = i;
	}
	bubbleSort(&sorted[1], dim);
	n_harm = max_harm = 0;
	if (max_harm > dim)
	    n_harm = dim;

	for(i=st; i<=st+dim; i++) {		
	    top = 0;
	    for(j=1; j<=n_harm; j++)
		if (model.A[i] == sorted[j].amp)
		    top = 1;
		
	    if (!top) {
		model.phi[i] = 0.0; /* make sure */
		//model.phi[i] = PI*(1.0 - 2.0*(float)rand()/RAND_MAX);
		model.phi[i] = phi_prev[i] + i*N*model.Wo;
		r++;
	    }
	    else
		nr++;
	}
	printf("r %d nr %d dim %d\n", r, nr, dim);
	    
	for(i=0; i<n_harm; i++)
	    ind[i] = sorted[i+1].index;
	for(i=n_harm; i<max_harm; i++)
	    ind[i] = 0;
	dump_hephase(ind, max_harm);
    }

    #ifdef UPPER
    for(i=3*model.L/4; i<=model.L; i++) {
	//model.phi[i] = phi_prev[i] + i*N*model.Wo;
	model.phi[i] = PI*(1.0 - 2.0*(float)rand()/RAND_MAX);
    }
    #endif
    for(i=1; i<=model.L; i++)
	phi_prev[i] = model.phi[i];	    
    #endif

    //#define HF
    #ifdef HF
    for(i=1; i<model.L/4; i++)
	model.phi[i] += (PI/8)*(1.0 - 2.0*(float)rand()/RAND_MAX);
	
    for(i=model.L/4; i<3*model.L/4; i++)
	model.phi[i] += (PI/4)*(1.0 - 2.0*(float)rand()/RAND_MAX);
	
    for(i=3*model.L/4; i<=model.L; i++)
	model.phi[i] += (PI/2)*(1.0 - 2.0*(float)rand()/RAND_MAX);	
    #endif

    //#define QUANT
    #ifdef QUANT
    for(i=1; i<=MAX_AMP; i++)
	model.phi[i] += (PI/4)*(1.0 - 2.0*(float)rand()/RAND_MAX);
    #endif

    //#define PRED
    #ifdef PRED
    for(i=1; i<=MAX_AMP; i++) {
	if ((frames % 2) != 0) {
	    /* predict on even frames */
	    model.phi[i] = phi_prev[i] + N*i*(model.Wo);
	}
	else {		
	    /* 2 bit quantise on odd */
	    model.phi[i] += (PI/8)*(1.0 - 2.0*(float)rand()/RAND_MAX);
	}
	phi_prev[i] = model.phi[i];	    
	Wo_prev = model.Wo;
    }	
    //if ((frames % 2) != 0)
    //    printf("frame: %d\n", frames);
    #endif

    //#define PRED_ERR
    #ifdef PRED_ERR
    for(i=model->L/4+1; i<=model->L/2; i++) {
	float pred = pexp->phi_prev[i] + N*i*(model->Wo);
	float err = pred - model->phi[i];
	err = atan2(sin(err),cos(err));
	printf("%f\n",err);
    }
    #endif

    //#define PRED2
    #ifdef PRED2
    /*
      if (fabs(model.Wo - Wo_prev) > 0.1*model.Wo) {
      for(i=1; i<=model.L/2; i++)
      phi_prev[i] = (PI/8)*(1.0 - 2.0*(float)rand()/RAND_MAX);
      }
      else
      printf("%d\n", frames);
    */
    #endif   

    #ifdef OLD
    if ((frames % 2) != 0) {
	/* predict on even frames */
	for(i=model.L/4+1; i<=model.L; i++)
	    model.phi[i] = phi_prev[i] + N*i*(model.Wo);
    }
    else {		
	/* 3 bit quantise on odd */
	for(i=model.L/4+1; i<=model.L; i++) {
	    float pred = phi_prev[i] + N*i*(model.Wo);
	    if (pred > model.phi[i]) 
		model.phi[i] = pred - PI/8;
	    else
		model.phi[i] = pred + PI/8;
		    
	    //model.phi[i] += (PI/8)*(1.0 - 2.0*(float)rand()/RAND_MAX);
	}
    }
    #endif

   #ifdef QUANT_SIM
   for(i=model.L/4+1; i<=model.L/2; i++) {
	float pred = phi_prev[i] + N*i*(model.Wo);
	float err = pred - model.phi[i];
	err = atan2(sin(err),cos(err));
	float interval = 0.25;
	int ind = floor(err/(interval*PI)+0.5);
	float qerr = ind*interval*PI;
	//printf("%f %d %f\n", err, ind, ind*0.25*PI);
	if (pred > model.phi[i]) 
	    model.phi[i] = pred - qerr;
	else
	    model.phi[i] = pred + qerr;
		    
    }

    for(i=model.L/2+1; i<=7*model.L/8; i++) {
	float pred = phi_prev[i] + N*i*(model.Wo);
	float err = pred - model.phi[i];
	err = atan2(sin(err),cos(err));
	float interval = 0.5;
	int ind = floor(err/(interval*PI)+0.5);
	float qerr = ind*interval*PI;
	//printf("%f %d %f\n", err, ind, ind*0.25*PI);
	if (pred > model.phi[i]) 
	    model.phi[i] = pred - qerr;
	else
	    model.phi[i] = pred + qerr;
		    
    }
    #endif

    /*
      for(i=1; i<=model.L/4; i++) {
      model.A[i] = 0.0;
      }
      for(i=model.L/2+1; i<=model.L; i++) {
      model.A[i] = 0.0;
      }
    */
	    

static void print_pred_error(struct PEXP *pexp, MODEL *model, int start, int end) {
    int i;
    for(i=start; i<=end; i++) {
	float pred = pexp->phi_prev[i] + N*i*(model->Wo);
	float err = pred - model->phi[i];
	err = atan2(sin(err),cos(err));
	printf("%f\n",err);
    }
}


static void predict_phases(struct PEXP *pexp, MODEL *model, int start, int end) {
    int i;

    for(i=start; i<=end; i++) {
	model->phi[i] = pexp->phi_prev[i] + N*i*(model->Wo + pexp->Wo_prev)/2.0;
	//model->phi[i] = pexp->phi_prev[i] + N*i*model->Wo;
    }
}

static void struct_phases(struct PEXP *pexp, MODEL *model, int start, int end) {
    int i;

    model->phi[1] = pexp->phi_prev[1] + N*model->Wo;

    for(i=start; i<=end; i++)
	model->phi[i] = model->phi[1]*i;
   
}

static void skip_phases(struct PEXP *pexp, MODEL *model, int start, int end) {
    int i;

    for(i=start; i<=end; i+=2)
	model->phi[i] = model->phi[i-1] - model->phi[i-2];
   
}

static void rand_phases(MODEL *model, int start, int end) {
    int i;

    for(i=start; i<=end; i++)
	model->phi[i] = PI*(1.0 - 2.0*(float)rand()/RAND_MAX);
   
}

static void quant_phase(float *phase, float min, float max, int bits) {
    int   levels = 1 << bits; 
    int   index;
    float norm, step;

    norm = (*phase - min)/(max - min);
    index = floor(levels*norm);

    //printf("phase %f norm %f index %d ", *phase, norm, index);
    if (index < 0 ) index = 0;
    if (index > (levels-1)) index = levels-1;
    //printf("index %d ", index);
    step = (max - min)/levels;
    *phase = min + step*index + 0.5*step;
    //printf("step %f phase %f\n", step, *phase);
}

static void quant_phases(MODEL *model, int start, int end, int bits) {
    int i;

    for(i=start; i<=end; i++) {
	quant_phase(&model->phi[i], -PI, PI, bits);
    }
}

static void fixed_bits_per_frame(struct PEXP *pexp, MODEL *model, int m, int budget) {
    int res, finished;

    res = 3;
    finished = 0;

    while(!finished) {
	if (m > model->L/2)
	    res = 2;
	if (((budget - res) < 0) || (m > model->L))
	    finished = 1;
	else {
	    quant_phase(&model->phi[m], -PI, PI, res);	    
	    budget -= res;
	    m++;
	}
    }
    printf("m: %d L: %d budget: %d\n", m, model->L, budget);
    predict_phases(pexp, model, m, model->L);
    //rand_phases(model, m, model->L);
}

/* used to plot histogram of quantisation error, for 3 bits, 8 levels,
   should be uniform between +/- PI/8 */

static void check_phase_quant(MODEL *model, float tol)
{
    int m;
    float phi_before[MAX_AMP];

    for(m=1; m<=model->L; m++)
	phi_before[m] = model->phi[m];

    quant_phases(model, 1, model->L, 3);

    for(m=1; m<=model->L; m++) {
	float err = phi_before[m] - model->phi[m];
	printf("%f\n", err);
	if (fabs(err) > tol)
	    exit(0);
    }
}


static void repeat_phases(MODEL *model, int period)
{
    int m,i;

    for(m=period+1,i=1; m<=model->L; m++,i++) {
	model->phi[m] = model->phi[m-1] + model->phi[i+1] - model->phi[i];
    }

}


static float est_phi1(MODEL *model, int start, int end)
{
    int m;
    float delta, s, c, phi1_est;

    if (end > model->L) 
	end = model->L;

    s = c = 0.0;
    for(m=start; m<end; m++) {
	delta = model->phi[m+1] - model->phi[m];
	s += model->A[m]*sin(delta);
	c += model->A[m]*cos(delta);
    }

    phi1_est = atan2(s,c);
    
    return phi1_est;
}

static void print_phi1_pred_error(MODEL *model, int start, int end)
{
    int m;
    float phi1_est;

    phi1_est = est_phi1(model, start, end);

    for(m=start; m<end; m++) {
	float err = model->phi[m+1] - model->phi[m] - phi1_est;
	//float err = model->phi[m] - phi1_est*m;
	err = atan2(sin(err),cos(err));
	printf("%f\n", err);
    }
}


static void first_order_band(MODEL *model, int start, int end, float phi1_est)
{
    int   m;
    float pred_err, av_pred_err;
    float c,s;

    s = c = 0.0;
    for(m=start; m<end; m++) {
	pred_err = model->phi[m] - phi1_est*m;
	s += model->A[m]*sin(pred_err);
	c += model->A[m]*cos(pred_err);
    }

    av_pred_err = atan2(s,c);
    for(m=start; m<end; m++) {
	model->phi[m] = av_pred_err + phi1_est*m;
	model->phi[m] = atan2(sin(model->phi[m]), cos(model->phi[m]));
    }

}


static void sub_linear(MODEL *model, int start, int end, float phi1_est)
{
    int   m;

    for(m=start; m<end; m++) {
	model->phi[m] = m*phi1_est;
    }
}


static void top_amp(struct PEXP *pexp, MODEL *model, int start, int end, int n_harm, int pred)
{
    int removed = 0, not_removed = 0;
    int top, i, j;
    struct AMPINDEX sorted[MAX_AMP];

    /* sort into acending order of amplitude */

    for(i=start,j=0; i<end; i++,j++) {
	sorted[j].amp = model->A[i];
	sorted[j].index = i;
    }
    bubbleSort(&sorted[1], end-start);

    /* keep phase of top n_harm, predict others */

    for(i=start; i<end; i++) {		
	top = 0;
	for(j=0; j<n_harm; j++)
	    if (model->A[i] == sorted[j].amp)
		top = 1;
		
	if (!top) {
	    model->phi[i] = 0.0; /* make sure */
	    if (pred)
		model->phi[i] = pexp->phi_prev[i] + i*N*(model->Wo + pexp->Wo_prev)/2.0;
	    else
		model->phi[i] = PI*(1.0 - 2.0*(float)rand()/RAND_MAX); // note: try rand for higher harms
	    removed++;
	}
	else {
	    /* need to make this work thru budget of bits */
	    //quant_phase(&model->phi[i], -PI, PI, 3);	    
	    not_removed++;
	}
    }
    //printf("dim: %d rem %d not_rem %d\n", end-start, removed, not_removed);
	    
}


/*---------------------------------------------------------------------------* \

  phase_experiment()

  Phase quantisation experiments.

\*---------------------------------------------------------------------------*/

void phase_experiment(struct PEXP *pexp, MODEL *model) {
    int m;
    float phi1_est;

    assert(pexp != NULL);

    //fixed_bits_per_frame(pexp, model, 40);
    //quant_phases(model, 1, model->L, 3);
    //print_pred_error(pexp, model, 1, model->L/2);
    //struct_phases(pexp, model, 1, model->L/4);
    //rand_phases(model, 10, model->L);
    //for(m=1; m<=model->L; m++)
    //	model->A[m] = 0.0;
    //model->A[model->L/2] = 1000;
    //repeat_phases(model, 20);
    //predict_phases(pexp, model, 1, model->L/4);
    //quant_phases(model, 1, 10, 3);
    //quant_phases(model, 10, 20, 2);
    //repeat_phases(model, 20);
    //rand_phases(model, 3*model->L/4, model->L);
    // print_phi1_pred_error(model, 1, model->L);
    //predict_phases(pexp, model, 1, model->L/4);
    //first_order_band(model, model->L/4, model->L/2);
    //first_order_band(model, model->L/2, 3*model->L/4);
    //if (fabs(model->Wo - pexp->Wo_prev)< 0.1*model->Wo)
    

    //phi1_est = est_phi1(model, 1, model->L/2);
    //first_order_band(model, 1, model->L/4, phi1_est);	
    //sub_linear(model, 1, model->L/4, phi1_est);

    //top_amp(pexp, model, 1, model->L/4, 4);
    //top_amp(pexp, model, model->L/4, model->L/2, 4);

    //first_order_band(model, 1, model->L/4, phi1_est);	
    //first_order_band(model, model->L/4, model->L/2, phi1_est);	

    //if (fabs(model->Wo - pexp->Wo_prev) > 0.2*model->Wo)
    //	rand_phases(model, model->L/2, model->L);
	
    //top_amp(pexp, model, 1, model->L/4, 4);
    //top_amp(pexp, model, model->L/4, model->L/2, 8);
    //top_amp(pexp, model, model->L/4+1, model->L/2, 10, 1);
    //top_amp(pexp, model, 1, model->L/4, 10, 1);
    //top_amp(pexp, model, model->L/4+1, 3*model->L/4, 10, 1);
    //top_amp(pexp, model, 1, 3*model->L/4, 20, 1);

    #ifdef REAS_CAND1
    predict_phases(pexp, model, 1, model->L/4);
    top_amp(pexp, model, model->L/4+1, 3*model->L/4, 10, 1);
    rand_phases(model, 3*model->L/4+1, model->L);
    #endif

    #ifdef REAS_CAND2
    if ((pexp->frames % 2) == 0) {
	//printf("quant\n");
	predict_phases(pexp, model, 1, model->L/4);	
	//top_amp(pexp, model, model->L/4+1, 3*model->L/4, 20, 1);
	top_amp(pexp, model,  model->L/4+1, 7*model->L/8, 20, 1);
	rand_phases(model, 7*model->L/8+1, model->L);
     }
    else {
	//printf("predict\n");
	predict_phases(pexp, model, 1, model->L);
    }
    #endif


    for(m=1; m<=model->L; m++)
    	model->phi[m] = atan2(sin(model->phi[m]),cos(model->phi[m]));

    /* update states */

    for(m=1; m<model->L; m++)
	pexp->phi_prev[m] = model->phi[m];	    
    pexp->Wo_prev = model->Wo;
    pexp->frames++;
}

