/*---------------------------------------------------------------------------*\
                                                                             
  FILE........: quantise.c
  AUTHOR......: David Rowe                                                     
  DATE CREATED: 31/5/92                                                       
                                                                             
  Quantisation functions for the sinusoidal coder.  
                                                                             
\*---------------------------------------------------------------------------*/

/*
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
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "defines.h"
#include "dump.h"
#include "quantise.h"
#include "lpc.h"
#include "lsp.h"
#include "fft.h"

#define LSP_DELTA1 0.01         /* grid spacing for LSP root searches */

/*---------------------------------------------------------------------------*\
									      
                          FUNCTION HEADERS

\*---------------------------------------------------------------------------*/

float speech_to_uq_lsps(float lsp[], float ak[], float Sn[], float w[], 
			int order);

/*---------------------------------------------------------------------------*\
									      
                             FUNCTIONS

\*---------------------------------------------------------------------------*/

int lsp_bits(int i) {
    return lsp_cb[i].log2m;
}

int lspd_bits(int i) {
    return lsp_cbd[i].log2m;
}

int lsp_pred_vq_bits(int i) {
    return lsp_cbjvm[i].log2m;
}

/*---------------------------------------------------------------------------*\

  quantise_init

  Loads the entire LSP quantiser comprised of several vector quantisers
  (codebooks).

\*---------------------------------------------------------------------------*/

void quantise_init()
{
}

/*---------------------------------------------------------------------------*\

  quantise

  Quantises vec by choosing the nearest vector in codebook cb, and
  returns the vector index.  The squared error of the quantised vector
  is added to se.

\*---------------------------------------------------------------------------*/

long quantise(const float * cb, float vec[], float w[], int k, int m, float *se)
/* float   cb[][K];	current VQ codebook		*/
/* float   vec[];	vector to quantise		*/
/* float   w[];         weighting vector                */
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
	    e += pow(diff*w[i],2.0);
	}
	if (e < beste) {
	    beste = e;
	    besti = j;
	}
   }

   *se += beste;

   return(besti);
}

/*---------------------------------------------------------------------------*\
									      
  lspd_quantise

  Scalar/VQ LSP difference quantiser.

\*---------------------------------------------------------------------------*/

void lspd_quantise(
  float lsp[], 
  float lsp_[],
  int   order
) 
{
    int   i,k,m,index;
    float lsp_hz[LPC_MAX];
    float lsp__hz[LPC_MAX];
    float dlsp[LPC_MAX];
    float dlsp_[LPC_MAX];
    float wt[LPC_MAX];
    const float *cb;
    float se;
    int   indexes[LPC_MAX];

    for(i=0; i<LPC_ORD; i++) {
	wt[i] = 1.0;
    }

    /* convert from radians to Hz so we can use human readable
       frequencies */

    for(i=0; i<order; i++)
	lsp_hz[i] = (4000.0/PI)*lsp[i];

    dlsp[0] = lsp_hz[0];
    for(i=1; i<order; i++)
    	dlsp[i] = lsp_hz[i] - lsp_hz[i-1];


    /* scalar quantisers for LSP differences 1..4 */

    wt[0] = 1.0;
    for(i=0; i<4; i++) {
	if (i) 
	    dlsp[i] = lsp_hz[i] - lsp__hz[i-1];	    
	else
	    dlsp[0] = lsp_hz[0];

	k = lsp_cbd[i].k;
	m = lsp_cbd[i].m;
	cb = lsp_cbd[i].cb;
	indexes[i] = quantise(cb, &dlsp[i], wt, k, m, &se);
 	dlsp_[i] = cb[indexes[i]*k];

	if (i) 
	    lsp__hz[i] = lsp__hz[i-1] + dlsp_[i];
	else
	    lsp__hz[0] = dlsp_[0];
	lsp_[i] = (PI/4000.0)*lsp__hz[i];
    }

    //#define PREV_VQ
#ifdef PREV_VQ
#define WGHT
#ifdef WGHT
    for(i=4; i<9; i++) {
	wt[i] = 1.0/(lsp[i]-lsp[i-1]) + 1.0/(lsp[i+1]-lsp[i]);
	//printf("wt[%d] = %f\n", i, wt[i]);
    }
    wt[9] = 1.0/(lsp[i]-lsp[i-1]);
#endif

    /* VQ LSPs 5,6,7,8,9,10 */

    ncb = 4;
    nlsp = 4;
    k = lsp_cbd[ncb].k;
    m = lsp_cbd[ncb].m;
    cb = lsp_cbd[ncb].cb;
    index = quantise(cb, &lsp[nlsp], &wt[nlsp], k, m, &se);
    lsp_[nlsp] = cb[index*k];
    lsp_[nlsp+1] = cb[index*k+1];
    lsp_[nlsp+2] = cb[index*k+2];
    lsp_[nlsp+3] = cb[index*k+3];
    lsp_[nlsp+4] = cb[index*k+4];
    lsp_[nlsp+5] = cb[index*k+5];

   /* convert back to radians */

    for(i=0; i<4; i++)
	lsp_[i] = (PI/4000.0)*lsp__hz[i];
#else
    /* VQ LSPs 5,6,7,8,9,10 */

    k = lsp_cbjnd[4].k;
    m = lsp_cbjnd[4].m;
    cb = lsp_cbjnd[4].cb;
    index = quantise(cb, &lsp_hz[4], &wt[4], k, m, &se);
    //printf("index = %4d: ", index);
    for(i=4; i<LPC_ORD; i++) {
	lsp_[i] = cb[index*k+i-4]*(PI/4000.0);
	//printf("%4.f (%4.f) ", lsp_hz[i], cb[index*k+i-4]);
    }
    //printf("\n");
#endif
}

/*---------------------------------------------------------------------------*\
									      
  lspvq_quantise

  Vector LSP quantiser.

\*---------------------------------------------------------------------------*/

void lspvq_quantise(
  float lsp[], 
  float lsp_[],
  int   order
) 
{
    int   i,k,m,ncb, nlsp;
    float  wt[LPC_ORD], lsp_hz[LPC_ORD];
    const float *cb;
    float se;
    int   index;

    for(i=0; i<LPC_ORD; i++) {
	wt[i] = 1.0;
	lsp_hz[i] = 4000.0*lsp[i]/PI;
    }

    /* scalar quantise LSPs 1,2,3,4 */

    /* simple uniform scalar quantisers */

   for(i=0; i<4; i++) {
	k = lsp_cb[i].k;
	m = lsp_cb[i].m;
	cb = lsp_cb[i].cb;
	index = quantise(cb, &lsp_hz[i], wt, k, m, &se);
	lsp_[i] = cb[index*k]*PI/4000.0;
    }

   //#define WGHT
#ifdef WGHT
    for(i=4; i<9; i++) {
	wt[i] = 1.0/(lsp[i]-lsp[i-1]) + 1.0/(lsp[i+1]-lsp[i]);
	//printf("wt[%d] = %f\n", i, wt[i]);
    }
    wt[9] = 1.0/(lsp[i]-lsp[i-1]);
#endif

    /* VQ LSPs 5,6,7,8,9,10 */

    ncb = 4;
    nlsp = 4;
    k = lsp_cbjnd[ncb].k;
    m = lsp_cbjnd[ncb].m;
    cb = lsp_cbjnd[ncb].cb;
    index = quantise(cb, &lsp_hz[nlsp], &wt[nlsp], k, m, &se);
    for(i=4; i<LPC_ORD; i++) {
	lsp_[i] = cb[index*k+i-4]*(PI/4000.0);
	//printf("%4.f (%4.f) ", lsp_hz[i], cb[index*k+i-4]);
    }
}

/*---------------------------------------------------------------------------*\
									      
  lspjnd_quantise

  Experimental JND LSP quantiser.

\*---------------------------------------------------------------------------*/

void lspjnd_quantise(float lsps[], float lsps_[], int order) 
{
    int   i,k,m;
    float  wt[LPC_ORD], lsps_hz[LPC_ORD];
    const float *cb;
    float se = 0.0;
    int   index;
 
    for(i=0; i<LPC_ORD; i++) {
	wt[i] = 1.0;
    }

    /* convert to Hz */

    for(i=0; i<LPC_ORD; i++) {
	lsps_hz[i] = lsps[i]*(4000.0/PI);
	lsps_[i] = lsps[i];
    }

    /* simple uniform scalar quantisers */

    for(i=0; i<4; i++) {
	k = lsp_cbjnd[i].k;
	m = lsp_cbjnd[i].m;
	cb = lsp_cbjnd[i].cb;
	index = quantise(cb, &lsps_hz[i], wt, k, m, &se);
	lsps_[i] = cb[index*k]*(PI/4000.0);
    }

    /* VQ LSPs 5,6,7,8,9,10 */

    k = lsp_cbjnd[4].k;
    m = lsp_cbjnd[4].m;
    cb = lsp_cbjnd[4].cb;
    index = quantise(cb, &lsps_hz[4], &wt[4], k, m, &se);
    //printf("k = %d m = %d c[0] %f cb[k] %f\n", k,m,cb[0],cb[k]);
    //printf("index = %4d: ", index);
    for(i=4; i<LPC_ORD; i++) {
	lsps_[i] = cb[index*k+i-4]*(PI/4000.0);
	//printf("%4.f (%4.f) ", lsps_hz[i], cb[index*k+i-4]);
    }
    //printf("\n");
}

void compute_weights(const float *x, float *w, int ndim);

/*---------------------------------------------------------------------------*\
									      
  lspdt_quantise

  LSP difference in time quantiser.  Split VQ, encoding LSPs 1-4 with
  one VQ, and LSPs 5-10 with a second.  Update of previous lsp memory
  is done outside of this function to handle dT between 10 or 20ms
  frames.

  mode        action
  ------------------

  LSPDT_ALL   VQ LSPs 1-4 and 5-10
  LSPDT_LOW   Just VQ LSPs 1-4, for LSPs 5-10 just copy previous
  LSPDT_HIGH  Just VQ LSPs 5-10, for LSPs 1-4 just copy previous

\*---------------------------------------------------------------------------*/

void lspdt_quantise(float lsps[], float lsps_[], float lsps__prev[], int mode)
{
    int   i;
    float wt[LPC_ORD];
    float lsps_dt[LPC_ORD];
#ifdef TRY_LSPDT_VQ
    int k,m;
    int   index;
    const float *cb;
    float se = 0.0;
#endif // TRY_LSPDT_VQ
    
    //compute_weights(lsps, wt, LPC_ORD);
    for(i=0; i<LPC_ORD; i++) {
    wt[i] = 1.0;
    }
 
    //compute_weights(lsps, wt, LPC_ORD );

    for(i=0; i<LPC_ORD; i++) {
	lsps_dt[i] = lsps[i] - lsps__prev[i];
	lsps_[i] = lsps__prev[i];
    }

    //#define TRY_LSPDT_VQ
#ifdef TRY_LSPDT_VQ
    /* this actually improves speech a bit, but 40ms updates works surprsingly well.... */
    k = lsp_cbdt[0].k;
    m = lsp_cbdt[0].m;
    cb = lsp_cbdt[0].cb;
    index = quantise(cb, lsps_dt, wt, k, m, &se);
    for(i=0; i<LPC_ORD; i++) {
	lsps_[i] += cb[index*k + i];
    }
#endif

}

#define MIN(a,b) ((a)<(b)?(a):(b))
#define MAX_ENTRIES 16384

void compute_weights(const float *x, float *w, int ndim)
{
  int i;
  w[0] = MIN(x[0], x[1]-x[0]);
  for (i=1;i<ndim-1;i++)
    w[i] = MIN(x[i]-x[i-1], x[i+1]-x[i]);
  w[ndim-1] = MIN(x[ndim-1]-x[ndim-2], PI-x[ndim-1]);
  
  for (i=0;i<ndim;i++)
    w[i] = 1./(.01+w[i]);
  //w[0]*=3;
  //w[1]*=2;
}

int find_nearest(const float *codebook, int nb_entries, float *x, int ndim)
{
  int i, j;
  float min_dist = 1e15;
  int nearest = 0;
  
  for (i=0;i<nb_entries;i++)
  {
    float dist=0;
    for (j=0;j<ndim;j++)
      dist += (x[j]-codebook[i*ndim+j])*(x[j]-codebook[i*ndim+j]);
    if (dist<min_dist)
    {
      min_dist = dist;
      nearest = i;
    }
  }
  return nearest;
}

int find_nearest_weighted(const float *codebook, int nb_entries, float *x, const float *w, int ndim)
{
  int i, j;
  float min_dist = 1e15;
  int nearest = 0;
  
  for (i=0;i<nb_entries;i++)
  {
    float dist=0;
    for (j=0;j<ndim;j++)
      dist += w[j]*(x[j]-codebook[i*ndim+j])*(x[j]-codebook[i*ndim+j]);
    if (dist<min_dist)
    {
      min_dist = dist;
      nearest = i;
    }
  }
  return nearest;
}

void lspjvm_quantise(float *x, float *xq, int ndim)
{
  int i, n1, n2, n3;
  float err[ndim], err2[ndim], err3[ndim];
  float w[ndim], w2[ndim], w3[ndim];
  const float *codebook1 = lsp_cbjvm[0].cb;
  const float *codebook2 = lsp_cbjvm[1].cb;
  const float *codebook3 = lsp_cbjvm[2].cb;

  w[0] = MIN(x[0], x[1]-x[0]);
  for (i=1;i<ndim-1;i++)
    w[i] = MIN(x[i]-x[i-1], x[i+1]-x[i]);
  w[ndim-1] = MIN(x[ndim-1]-x[ndim-2], PI-x[ndim-1]);
  
  compute_weights(x, w, ndim);
  
  n1 = find_nearest(codebook1, lsp_cbjvm[0].m, x, ndim);
  
  for (i=0;i<ndim;i++)
  {
    xq[i] = codebook1[ndim*n1+i];
    err[i] = x[i] - xq[i];
  }
  for (i=0;i<ndim/2;i++)
  {
    err2[i] = err[2*i];  
    err3[i] = err[2*i+1];
    w2[i] = w[2*i];  
    w3[i] = w[2*i+1];
  }
  n2 = find_nearest_weighted(codebook2, lsp_cbjvm[1].m, err2, w2, ndim/2);
  n3 = find_nearest_weighted(codebook3, lsp_cbjvm[2].m, err3, w3, ndim/2);
  
  for (i=0;i<ndim/2;i++)
  {
    xq[2*i] += codebook2[ndim*n2/2+i];
    xq[2*i+1] += codebook3[ndim*n3/2+i];
  }
}

void check_lsp_order(float lsp[], int lpc_order)
{
    int   i;
    float tmp;

    for(i=1; i<lpc_order; i++)
	if (lsp[i] < lsp[i-1]) {
	    //printf("swap %d\n",i);
	    tmp = lsp[i-1];
	    lsp[i-1] = lsp[i]-0.05;
	    lsp[i] = tmp+0.05;
	}
}

void force_min_lsp_dist(float lsp[], int lpc_order)
{
    int   i;

    for(i=1; i<lpc_order; i++)
	if ((lsp[i]-lsp[i-1]) < 0.01) {
	    lsp[i] += 0.01;
	}
}

/*---------------------------------------------------------------------------*\
									      
  lpc_model_amplitudes

  Derive a LPC model for amplitude samples then estimate amplitude samples
  from this model with optional LSP quantisation.

  Returns the spectral distortion for this frame.

\*---------------------------------------------------------------------------*/

float lpc_model_amplitudes(
  float  Sn[],			/* Input frame of speech samples */
  float  w[],			
  MODEL *model,			/* sinusoidal model parameters */
  int    order,                 /* LPC model order */
  int    lsp_quant,             /* optional LSP quantisation if non-zero */
  float  ak[]                   /* output aks */
)
{
  float Wn[M];
  float R[LPC_MAX+1];
  float E;
  int   i,j;
  float snr;	
  float lsp[LPC_MAX];
  float lsp_hz[LPC_MAX];
  float lsp_[LPC_MAX];
  int   roots;                  /* number of LSP roots found */
  float wt[LPC_MAX];

  for(i=0; i<M; i++)
    Wn[i] = Sn[i]*w[i];
  autocorrelate(Wn,R,M,order);
  levinson_durbin(R,ak,order);
  
  E = 0.0;
  for(i=0; i<=order; i++)
      E += ak[i]*R[i];
 
  for(i=0; i<order; i++)
      wt[i] = 1.0;

  if (lsp_quant) {
    roots = lpc_to_lsp(ak, order, lsp, 5, LSP_DELTA1);
    if (roots != order)
	printf("LSP roots not found\n");

    /* convert from radians to Hz to make quantisers more
       human readable */

    for(i=0; i<order; i++)
	lsp_hz[i] = (4000.0/PI)*lsp[i];
    
#ifdef NOT_NOW
    /* simple uniform scalar quantisers */

    for(i=0; i<10; i++) {
	k = lsp_cb[i].k;
	m = lsp_cb[i].m;
	cb = lsp_cb[i].cb;
	index = quantise(cb, &lsp_hz[i], wt, k, m, &se);
	lsp_hz[i] = cb[index*k];
    }
#endif
    
    /* experiment: simulating uniform quantisation error
    for(i=0; i<order; i++)
	lsp[i] += PI*(12.5/4000.0)*(1.0 - 2.0*(float)rand()/RAND_MAX);
    */

    for(i=0; i<order; i++)
	lsp[i] = (PI/4000.0)*lsp_hz[i];

    /* Bandwidth Expansion (BW).  Prevents any two LSPs getting too
       close together after quantisation.  We know from experiment
       that LSP quantisation errors < 12.5Hz (25Hz step size) are
       inaudible so we use that as the minimum LSP separation.
    */

    for(i=1; i<5; i++) {
	if (lsp[i] - lsp[i-1] < PI*(12.5/4000.0))
	    lsp[i] = lsp[i-1] + PI*(12.5/4000.0);
    }

    /* as quantiser gaps increased, larger BW expansion was required
       to prevent twinkly noises */

    for(i=5; i<8; i++) {
	if (lsp[i] - lsp[i-1] < PI*(25.0/4000.0))
	    lsp[i] = lsp[i-1] + PI*(25.0/4000.0);
    }
    for(i=8; i<order; i++) {
	if (lsp[i] - lsp[i-1] < PI*(75.0/4000.0))
	    lsp[i] = lsp[i-1] + PI*(75.0/4000.0);
    }

    for(j=0; j<order; j++) 
	lsp_[j] = lsp[j];

    lsp_to_lpc(lsp_, ak, order);
#ifdef DUMP
    dump_lsp(lsp);
#endif
  }

#ifdef DUMP
  dump_E(E);
#endif
  #ifdef SIM_QUANT
  /* simulated LPC energy quantisation */
  {
      float e = 10.0*log10(E);
      e += 2.0*(1.0 - 2.0*(float)rand()/RAND_MAX);
      E = pow(10.0,e/10.0);
  }
  #endif

  aks_to_M2(ak,order,model,E,&snr, 1);   /* {ak} -> {Am} LPC decode */

  return snr;
}

/*---------------------------------------------------------------------------*\
                                                                         
   aks_to_M2()                                                             
                                                                         
   Transforms the linear prediction coefficients to spectral amplitude    
   samples.  This function determines A(m) from the average energy per    
   band using an FFT.                                                     
                                                                        
\*---------------------------------------------------------------------------*/

void aks_to_M2(
  float  ak[],	/* LPC's */
  int    order,
  MODEL *model,	/* sinusoidal model parameters for this frame */
  float  E,	/* energy term */
  float *snr,	/* signal to noise ratio for this frame in dB */
  int    dump   /* true to dump sample to dump file */
)
{
  COMP Pw[FFT_DEC];	/* power spectrum */
  int i,m;		/* loop variables */
  int am,bm;		/* limits of current band */
  float r;		/* no. rads/bin */
  float Em;		/* energy in band */
  float Am;		/* spectral amplitude sample */
  float signal, noise;

  r = TWO_PI/(FFT_DEC);

  /* Determine DFT of A(exp(jw)) --------------------------------------------*/

  for(i=0; i<FFT_DEC; i++) {
    Pw[i].real = 0.0;
    Pw[i].imag = 0.0; 
  }

  for(i=0; i<=order; i++)
    Pw[i].real = ak[i];
  fft(&Pw[0].real,FFT_DEC,1);

  /* Determine power spectrum P(w) = E/(A(exp(jw))^2 ------------------------*/

  for(i=0; i<FFT_DEC/2; i++)
    Pw[i].real = E/(Pw[i].real*Pw[i].real + Pw[i].imag*Pw[i].imag);
#ifdef DUMP
  if (dump) 
      dump_Pw(Pw);
#endif

  /* Determine magnitudes by linear interpolation of P(w) -------------------*/

  signal = noise = 0.0;
  for(m=1; m<=model->L; m++) {
    am = floor((m - 0.5)*model->Wo/r + 0.5);
    bm = floor((m + 0.5)*model->Wo/r + 0.5);
    Em = 0.0;

    for(i=am; i<bm; i++)
      Em += Pw[i].real;
    Am = sqrt(Em);

    signal += pow(model->A[m],2.0);
    noise  += pow(model->A[m] - Am,2.0);
    model->A[m] = Am;
  }
  *snr = 10.0*log10(signal/noise);
}

/*---------------------------------------------------------------------------*\
                                                       
  FUNCTION....: encode_Wo()	     
  AUTHOR......: David Rowe			      
  DATE CREATED: 22/8/2010 

  Encodes Wo using a WO_LEVELS quantiser.

\*---------------------------------------------------------------------------*/

int encode_Wo(float Wo)
{
    int   index;
    float Wo_min = TWO_PI/P_MAX;
    float Wo_max = TWO_PI/P_MIN;
    float norm;

    norm = (Wo - Wo_min)/(Wo_max - Wo_min);
    index = floor(WO_LEVELS * norm + 0.5);
    if (index < 0 ) index = 0;
    if (index > (WO_LEVELS-1)) index = WO_LEVELS-1;

    return index;
}

/*---------------------------------------------------------------------------*\
                                                       
  FUNCTION....: decode_Wo()	     
  AUTHOR......: David Rowe			      
  DATE CREATED: 22/8/2010 

  Decodes Wo using a WO_LEVELS quantiser.

\*---------------------------------------------------------------------------*/

float decode_Wo(int index)
{
    float Wo_min = TWO_PI/P_MAX;
    float Wo_max = TWO_PI/P_MIN;
    float step;
    float Wo;

    step = (Wo_max - Wo_min)/WO_LEVELS;
    Wo   = Wo_min + step*(index);

    return Wo;
}

/*---------------------------------------------------------------------------*\
                                                       
  FUNCTION....: encode_Wo_dt()	     
  AUTHOR......: David Rowe			      
  DATE CREATED: 6 Nov 2011 

  Encodes Wo difference from last frame.

\*---------------------------------------------------------------------------*/

int encode_Wo_dt(float Wo, float prev_Wo)
{
    int   index, mask, max_index, min_index;
    float Wo_min = TWO_PI/P_MAX;
    float Wo_max = TWO_PI/P_MIN;
    float norm;

    norm = (Wo - prev_Wo)/(Wo_max - Wo_min);
    index = floor(WO_LEVELS * norm + 0.5);
    //printf("ENC index: %d ", index);

    /* hard limit */
    
    max_index = (1 << (WO_DT_BITS-1)) - 1;
    min_index = - (max_index+1);
    if (index > max_index) index = max_index;
    if (index < min_index) index = min_index;
    //printf("max_index: %d  min_index: %d hard index: %d ",
    //	   max_index,  min_index, index);

    /* mask so that only LSB WO_DT_BITS remain, bit WO_DT_BITS is the sign bit */

    mask = ((1 << WO_DT_BITS) - 1);
    index &= mask;
    //printf("mask: 0x%x index: 0x%x\n", mask, index);

    return index;
}

/*---------------------------------------------------------------------------*\
                                                       
  FUNCTION....: decode_Wo_dt()	     
  AUTHOR......: David Rowe			      
  DATE CREATED: 6 Nov 2011 

  Decodes Wo using WO_DT_BITS difference from last frame.

\*---------------------------------------------------------------------------*/

float decode_Wo_dt(int index, float prev_Wo)
{
    float Wo_min = TWO_PI/P_MAX;
    float Wo_max = TWO_PI/P_MIN;
    float step;
    float Wo;
    int   mask;

    /* sign extend index */
    
    //printf("DEC index: %d ");
    if (index & (1 << (WO_DT_BITS-1))) {
	mask = ~((1 << WO_DT_BITS) - 1);
	index |= mask;
    }
    //printf("DEC mask: 0x%x  index: %d \n", mask, index);
    
    step = (Wo_max - Wo_min)/WO_LEVELS;
    Wo   = prev_Wo + step*(index);

    /* bit errors can make us go out of range leading to all sorts of
       probs like seg faults */

    if (Wo > Wo_max) Wo = Wo_max;
    if (Wo < Wo_min) Wo = Wo_min;

    return Wo;
}

/*---------------------------------------------------------------------------*\
                                                       
  FUNCTION....: speech_to_uq_lsps()	     
  AUTHOR......: David Rowe			      
  DATE CREATED: 22/8/2010 

  Analyse a windowed frame of time domain speech to determine LPCs
  which are the converted to LSPs for quantisation and transmission
  over the channel.

\*---------------------------------------------------------------------------*/

float speech_to_uq_lsps(float lsp[],
			float ak[],
		        float Sn[], 
		        float w[],
		        int   order
)
{
    int   i, roots;
    float Wn[M];
    float R[LPC_MAX+1];
    float E;

    for(i=0; i<M; i++)
	Wn[i] = Sn[i]*w[i];
    autocorrelate(Wn, R, M, order);
    levinson_durbin(R, ak, order);
  
    E = 0.0;
    for(i=0; i<=order; i++)
	E += ak[i]*R[i];
 
    /* 15 Hz BW expansion as I can't hear the difference and it may help
       help occasional fails in the LSP root finding.  Important to do this
       after energy caluclation to avoid -ve energy values.
    */

    for(i=0; i<=order; i++)
	ak[i] *= pow(0.994,(float)i);

    roots = lpc_to_lsp(ak, order, lsp, 5, LSP_DELTA1);
    if (roots != order) {
	/* use some benign LSP values we can use instead */
	for(i=0; i<order; i++)
	    lsp[i] = (PI/order)*(float)i;
    }

    return E;
}

/*---------------------------------------------------------------------------*\
                                                       
  FUNCTION....: encode_lsps_scalar()	     
  AUTHOR......: David Rowe			      
  DATE CREATED: 22/8/2010 

  Thirty-six bit sclar LSP quantiser. From a vector of unquantised
  (floating point) LSPs finds the quantised LSP indexes.

\*---------------------------------------------------------------------------*/

void encode_lsps_scalar(int indexes[], float lsp[], int order)
{
    int    i,k,m;
    float  wt[1];
    float  lsp_hz[LPC_MAX];
    const float * cb;
    float se;

    /* convert from radians to Hz so we can use human readable
       frequencies */

    for(i=0; i<order; i++)
	lsp_hz[i] = (4000.0/PI)*lsp[i];
    
    /* simple uniform scalar quantisers */

    wt[0] = 1.0;
    for(i=0; i<order; i++) {
	k = lsp_cb[i].k;
	m = lsp_cb[i].m;
	cb = lsp_cb[i].cb;
	indexes[i] = quantise(cb, &lsp_hz[i], wt, k, m, &se);
    }
}

/*---------------------------------------------------------------------------*\
                                                       
  FUNCTION....: decode_lsps_scalar()	     
  AUTHOR......: David Rowe			      
  DATE CREATED: 22/8/2010 

  From a vector of quantised LSP indexes, returns the quantised
  (floating point) LSPs.

\*---------------------------------------------------------------------------*/

void decode_lsps_scalar(float lsp[], int indexes[], int order)
{
    int    i,k;
    float  lsp_hz[LPC_MAX];
    const float * cb;

    for(i=0; i<order; i++) {
	k = lsp_cb[i].k;
	cb = lsp_cb[i].cb;
	lsp_hz[i] = cb[indexes[i]*k];
    }

    /* convert back to radians */

    for(i=0; i<order; i++)
	lsp[i] = (PI/4000.0)*lsp_hz[i];
}


/*---------------------------------------------------------------------------*\
                                                       
  FUNCTION....: encode_lsps_diff_freq_vq()	     
  AUTHOR......: David Rowe			      
  DATE CREATED: 15 November 2011

  Twenty-five bit LSP quantiser.  LSPs 1-4 are quantised with scalar
  LSP differences (in frequency, i.e difference from the previous
  LSP).  LSPs 5-10 are quantised with a VQ trained generated using
  vqtrainjnd.c

\*---------------------------------------------------------------------------*/

void encode_lsps_diff_freq_vq(int indexes[], float lsp[], int order)
{
    int    i,k,m;
    float  lsp_hz[LPC_MAX];
    float lsp__hz[LPC_MAX];
    float dlsp[LPC_MAX];
    float dlsp_[LPC_MAX];
    float wt[LPC_MAX];
    const float * cb;
    float se;

    for(i=0; i<LPC_ORD; i++) {
	wt[i] = 1.0;
    }

    /* convert from radians to Hz so we can use human readable
       frequencies */

    for(i=0; i<order; i++)
	lsp_hz[i] = (4000.0/PI)*lsp[i];
    
    /* scalar quantisers for LSP differences 1..4 */

    wt[0] = 1.0;
    for(i=0; i<4; i++) {
	if (i) 
	    dlsp[i] = lsp_hz[i] - lsp__hz[i-1];	    
	else
	    dlsp[0] = lsp_hz[0];

	k = lsp_cbd[i].k;
	m = lsp_cbd[i].m;
	cb = lsp_cbd[i].cb;
	indexes[i] = quantise(cb, &dlsp[i], wt, k, m, &se);
 	dlsp_[i] = cb[indexes[i]*k];

	if (i) 
	    lsp__hz[i] = lsp__hz[i-1] + dlsp_[i];
	else
	    lsp__hz[0] = dlsp_[0];
    }

    /* VQ LSPs 5,6,7,8,9,10 */

    k = lsp_cbjnd[4].k;
    m = lsp_cbjnd[4].m;
    cb = lsp_cbjnd[4].cb;
    indexes[4] = quantise(cb, &lsp_hz[4], &wt[4], k, m, &se);
}


/*---------------------------------------------------------------------------*\
                                                       
  FUNCTION....: decode_lsps_diff_freq_vq()	     
  AUTHOR......: David Rowe			      
  DATE CREATED: 15 Nov 2011

  From a vector of quantised LSP indexes, returns the quantised
  (floating point) LSPs.

\*---------------------------------------------------------------------------*/

void decode_lsps_diff_freq_vq(float lsp_[], int indexes[], int order)
{
    int    i,k,m;
    float  dlsp_[LPC_MAX];
    float  lsp__hz[LPC_MAX];
    const float * cb;

    /* scalar LSP differences */

    for(i=0; i<4; i++) {
	cb = lsp_cbd[i].cb;
	dlsp_[i] = cb[indexes[i]];
	if (i) 
	    lsp__hz[i] = lsp__hz[i-1] + dlsp_[i];
	else
	    lsp__hz[0] = dlsp_[0];
    }

    /* VQ */

    k = lsp_cbjnd[4].k;
    m = lsp_cbjnd[4].m;
    cb = lsp_cbjnd[4].cb;
    for(i=4; i<order; i++)
	lsp__hz[i] = cb[indexes[4]*k+i-4];

    /* convert back to radians */

    for(i=0; i<order; i++)
	lsp_[i] = (PI/4000.0)*lsp__hz[i];
}


/*---------------------------------------------------------------------------*\
                                                       
  FUNCTION....: encode_lsps_diff_time_vq()	     
  AUTHOR......: David Rowe			      
  DATE CREATED: 15 November 2011

  Seven bit LSP quantiser.  The difference in LSPs 1-4 from a previous
  frame are quantised with a hand constructed Just Noticable Differnce
  (JND) VQ.  LSPs 5-10 are just repeated from the last frame.

\*---------------------------------------------------------------------------*/

void encode_lsps_diff_time_vq(int indexes[], 
			      float lsps[], 
			      float lsps__prev[], 
			      int order)
{
    int    i,k,m;
    float  lsps_dt[LPC_ORD];
    float  wt[LPC_MAX];
    const  float * cb;
    float  se;

    for(i=0; i<LPC_ORD; i++) {
	wt[i] = 1.0;
    }

    /* Determine difference in time and convert from radians to Hz so
       we can use human readable frequencies */

    for(i=0; i<LPC_ORD; i++) {
	lsps_dt[i] = (4000/PI)*(lsps[i] - lsps__prev[i]);
    }
    
    /* VQ LSP dTs 1 to 4 */

    k = lsp_cbdt[0].k;
    m = lsp_cbdt[0].m;
    cb = lsp_cbdt[0].cb;
    *indexes = quantise(cb, lsps_dt, wt, k, m, &se);
}


/*---------------------------------------------------------------------------*\
                                                       
  FUNCTION....: decode_lsps_diff_time_vq()	     
  AUTHOR......: David Rowe			      
  DATE CREATED: 15 Nov 2011

  From a quantised LSP indexes, returns the quantised
  (floating point) LSPs.

\*---------------------------------------------------------------------------*/

void decode_lsps_diff_time_vq(
			      float lsps_[], 
			      int indexes[], 
			      float lsps__prev[],
			      int order)
{
    int    i,k,m;
    const  float * cb;

    for(i=0; i<order; i++)
	lsps_[i] = lsps__prev[i];

    k = lsp_cbdt[0].k;
    m = lsp_cbdt[0].m;
    cb = lsp_cbdt[0].cb;

    for(i=0; i<4; i++) {
	lsps_[i] += (PI/4000.0)*cb[*indexes*k + i];
    }

}


/*---------------------------------------------------------------------------*\
                                                       
  FUNCTION....: encode_lsps_vq()	     
  AUTHOR......: David Rowe			      
  DATE CREATED: 15 Feb 2012

  Multi-stage VQ LSP quantiser developed by Jean-Marc Valin.

\*---------------------------------------------------------------------------*/

void encode_lsps_vq(int *indexes, float *x, float *xq, int ndim)
{
  int i, n1, n2, n3;
  float err[ndim], err2[ndim], err3[ndim];
  float w[ndim], w2[ndim], w3[ndim];
  const float *codebook1 = lsp_cbjvm[0].cb;
  const float *codebook2 = lsp_cbjvm[1].cb;
  const float *codebook3 = lsp_cbjvm[2].cb;

  w[0] = MIN(x[0], x[1]-x[0]);
  for (i=1;i<ndim-1;i++)
    w[i] = MIN(x[i]-x[i-1], x[i+1]-x[i]);
  w[ndim-1] = MIN(x[ndim-1]-x[ndim-2], PI-x[ndim-1]);
  
  compute_weights(x, w, ndim);
  
  n1 = find_nearest(codebook1, lsp_cbjvm[0].m, x, ndim);
  
  for (i=0;i<ndim;i++)
  {
    xq[i]  = codebook1[ndim*n1+i];
    err[i] = x[i] - xq[i];
  }
  for (i=0;i<ndim/2;i++)
  {
    err2[i] = err[2*i];  
    err3[i] = err[2*i+1];
    w2[i] = w[2*i];  
    w3[i] = w[2*i+1];
  }
  n2 = find_nearest_weighted(codebook2, lsp_cbjvm[1].m, err2, w2, ndim/2);
  n3 = find_nearest_weighted(codebook3, lsp_cbjvm[2].m, err3, w3, ndim/2);
  
  indexes[0] = n1;
  indexes[1] = n2;
  indexes[2] = n3;
}


/*---------------------------------------------------------------------------*\
                                                       
  FUNCTION....: decode_lsps_vq()	     
  AUTHOR......: David Rowe			      
  DATE CREATED: 15 Feb 2012

\*---------------------------------------------------------------------------*/

void decode_lsps_vq(int *indexes, float *xq, int ndim)
{
  int i, n1, n2, n3;
  const float *codebook1 = lsp_cbjvm[0].cb;
  const float *codebook2 = lsp_cbjvm[1].cb;
  const float *codebook3 = lsp_cbjvm[2].cb;

  n1 = indexes[0];
  n2 = indexes[1];
  n3 = indexes[2];

  for (i=0;i<ndim;i++)
  {
    xq[i] = codebook1[ndim*n1+i];
  }
  for (i=0;i<ndim/2;i++)
  {
    xq[2*i] += codebook2[ndim*n2/2+i];
    xq[2*i+1] += codebook3[ndim*n3/2+i];
  }
}


/*---------------------------------------------------------------------------*\
                                                       
  FUNCTION....: bw_expand_lsps()	     
  AUTHOR......: David Rowe			      
  DATE CREATED: 22/8/2010 

  Applies Bandwidth Expansion (BW) to a vector of LSPs.  Prevents any
  two LSPs getting too close together after quantisation.  We know
  from experiment that LSP quantisation errors < 12.5Hz (25Hz step
  size) are inaudible so we use that as the minimum LSP separation.

\*---------------------------------------------------------------------------*/

void bw_expand_lsps(float lsp[],
		    int   order
)
{
    int i;

    for(i=1; i<4; i++) {
	
	if ((lsp[i] - lsp[i-1]) < 50*(PI/4000.0))
	    lsp[i] = lsp[i-1] + 50.0*(PI/4000.0);
	
    }

    /* As quantiser gaps increased, larger BW expansion was required
       to prevent twinkly noises.  This may need more experiment for
       different quanstisers.
    */

    for(i=4; i<order; i++) {
	if (lsp[i] - lsp[i-1] < PI*(100.0/4000.0))
	    lsp[i] = lsp[i-1] + PI*(100.0/4000.0);
    }
}

/*---------------------------------------------------------------------------*\
                                                       
  FUNCTION....: locate_lsps_jnd_steps()	     
  AUTHOR......: David Rowe			      
  DATE CREATED: 27/10/2011 

  Applies a form of Bandwidth Expansion (BW) to a vector of LSPs.
  Listening tests have determined that "quantising" the position of
  each LSP to the non-linear steps below introduces a "just noticable
  difference" in the synthesised speech.

  This operation can be used before quantisation to limit the input
  data to the quantiser to a number of discrete steps.

  This operation can also be used during quantisation as a form of
  hysteresis in the calculation of quantiser error.  For example if
  the quantiser target of lsp1 is 500 Hz, candidate vectors with lsp1
  of 515 and 495 Hz sound effectively the same.

\*---------------------------------------------------------------------------*/

void locate_lsps_jnd_steps(float lsps[], int order)
{
    int   i;
    float lsp_hz, step;

    assert(order == 10);

    /* quantise to 25Hz steps */
	    
    step = 25;
    for(i=0; i<2; i++) {
	lsp_hz = lsps[i]*4000.0/PI;
	lsp_hz = floor(lsp_hz/step + 0.5)*step;
	lsps[i] = lsp_hz*PI/4000.0;
	if (i) {
	    if (lsps[i] == lsps[i-1])
		lsps[i]   += step*PI/4000.0;

	}
    }

    /* quantise to 50Hz steps */

    step = 50;
    for(i=2; i<4; i++) {
	lsp_hz = lsps[i]*4000.0/PI;
	lsp_hz = floor(lsp_hz/step + 0.5)*step;
	lsps[i] = lsp_hz*PI/4000.0;
	if (i) {
	    if (lsps[i] == lsps[i-1])
		lsps[i] += step*PI/4000.0;

	}
    }

    /* quantise to 100Hz steps */

    step = 100;
    for(i=4; i<10; i++) {
	lsp_hz = lsps[i]*4000.0/PI;
	lsp_hz = floor(lsp_hz/step + 0.5)*step;
	lsps[i] = lsp_hz*PI/4000.0;
	if (i) {
	    if (lsps[i] == lsps[i-1])
		lsps[i] += step*PI/4000.0;

	}
    }
}


/*---------------------------------------------------------------------------*\
                                                       
  FUNCTION....: apply_lpc_correction()	     
  AUTHOR......: David Rowe			      
  DATE CREATED: 22/8/2010 

  Apply first harmonic LPC correction at decoder.  This helps improve
  low pitch males after LPC modelling, like hts1a and morig.

\*---------------------------------------------------------------------------*/

void apply_lpc_correction(MODEL *model)
{
    if (model->Wo < (PI*150.0/4000)) {
	model->A[1] *= 0.032;
    }
}

/*---------------------------------------------------------------------------*\
                                                       
  FUNCTION....: encode_energy()	     
  AUTHOR......: David Rowe			      
  DATE CREATED: 22/8/2010 

  Encodes LPC energy using an E_LEVELS quantiser.

\*---------------------------------------------------------------------------*/

int encode_energy(float e)
{
    int   index;
    float e_min = E_MIN_DB;
    float e_max = E_MAX_DB;
    float norm;

    e = 10.0*log10(e);
    norm = (e - e_min)/(e_max - e_min);
    index = floor(E_LEVELS * norm + 0.5);
    if (index < 0 ) index = 0;
    if (index > (E_LEVELS-1)) index = E_LEVELS-1;

    return index;
}

/*---------------------------------------------------------------------------*\
                                                       
  FUNCTION....: decode_energy()	     
  AUTHOR......: David Rowe			      
  DATE CREATED: 22/8/2010 

  Decodes energy using a E_LEVELS quantiser.

\*---------------------------------------------------------------------------*/

float decode_energy(int index)
{
    float e_min = E_MIN_DB;
    float e_max = E_MAX_DB;
    float step;
    float e;

    step = (e_max - e_min)/E_LEVELS;
    e    = e_min + step*(index);
    e    = pow(10.0,e/10.0);

    return e;
}

/*---------------------------------------------------------------------------*\
                                                       
  FUNCTION....: decode_amplitudes()	     
  AUTHOR......: David Rowe			      
  DATE CREATED: 22/8/2010 

  Given the amplitude quantiser indexes recovers the harmonic
  amplitudes.

\*---------------------------------------------------------------------------*/

float decode_amplitudes(MODEL *model, 
			float  ak[],
		        int    lsp_indexes[], 
		        int    energy_index,
			float  lsps[],
			float *e
)
{
    float snr;

    decode_lsps_scalar(lsps, lsp_indexes, LPC_ORD);
    bw_expand_lsps(lsps, LPC_ORD);
    lsp_to_lpc(lsps, ak, LPC_ORD);
    *e = decode_energy(energy_index);
    aks_to_M2(ak, LPC_ORD, model, *e, &snr, 1); 
    apply_lpc_correction(model);

    return snr;
}


static float ge_coeff[2] = {0.8, 0.9};

void compute_weights2(const float *x, const float *xp, float *w, int ndim)
{
  w[0] = 30;
  w[1] = 1;
  if (x[1]<0)
  {
     w[0] *= .6;
     w[1] *= .3;
  }
  if (x[1]<-10)
  {
     w[0] *= .3;
     w[1] *= .3;
  }
  /* Higher weight if pitch is stable */
  if (fabs(x[0]-xp[0])<.2)
  {
     w[0] *= 2;
     w[1] *= 1.5;
  } else if (fabs(x[0]-xp[0])>.5) /* Lower if not stable */
  {
     w[0] *= .5;
  }

  /* Lower weight for low energy */
  if (x[1] < xp[1]-10)
  {
     w[1] *= .5;
  }
  if (x[1] < xp[1]-20)
  {
     w[1] *= .5;
  }

  //w[0] = 30;
  //w[1] = 1;
  
  /* Square the weights because it's applied on the squared error */
  w[0] *= w[0];
  w[1] *= w[1];

}

/*---------------------------------------------------------------------------*\
                                                       
  FUNCTION....: quantize_ge()	     
  AUTHOR......: Jean-Marc Valin & David Rowe			      
  DATE CREATED: 29 Feb 2012

  Experimental joint Wo and LPC energy vector quantiser developed my
  Jean-Marc Valin.
  
\*---------------------------------------------------------------------------*/

void ge_quantise(float *x, float *xq)
{
  int          i, n1;
  float        err[2];
  float        w[2];
  const float *codebook1 = ge_cb[0].cb;
  int          nb_entries = ge_cb[0].m;
  int          ndim = ge_cb[0].k;

  //printf("ndim %d nb_entries %d\n", ndim, nb_entries);
  compute_weights2(x, xq, w, ndim);
  //exit(0);
  for (i=0;i<ndim;i++)
    err[i] = x[i]-ge_coeff[i]*xq[i];
  n1 = find_nearest_weighted(codebook1, nb_entries, err, w, ndim);
  
  for (i=0;i<ndim;i++)
  {
    xq[i] = ge_coeff[i]*xq[i] + codebook1[ndim*n1+i];
    err[i] -= codebook1[ndim*n1+i];
  }

}

