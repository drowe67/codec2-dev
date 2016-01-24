/*---------------------------------------------------------------------------*\

  FILE........: fsk.c
  AUTHOR......: Brady O'Brien
  DATE CREATED: 7 January 2016

  C Implementation of 2FSK modulator/demodulator, based on octave/fsk_horus.m

\*---------------------------------------------------------------------------*/

/*
  Copyright (C) 2016 David Rowe

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

/*---------------------------------------------------------------------------*\

                               DEFINES

\*---------------------------------------------------------------------------*/

/* P oversampling rate constant -- should probably be init-time configurable */
#define ct_P 8

/* Define this to enable EbNodB estimate */
/* This needs square roots, may take more cpu time than it's worth */
#define EST_EBNO

/*---------------------------------------------------------------------------*\

                               INCLUDES

\*---------------------------------------------------------------------------*/

#include <assert.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>

#include "fsk.h"
#include "comp_prim.h"
#include "kiss_fftr.h"
#include "modem_probe.h"

/*---------------------------------------------------------------------------*\

                               FUNCTIONS

\*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*\

  FUNCTION....: fsk_create
  AUTHOR......: Brady O'Brien
  DATE CREATED: 7 January 2016
  
  Create and initialize an instance of the FSK modem. Returns a pointer
  to the modem state/config struct. One modem config struct may be used
  for both mod and demod. returns NULL on failure.

\*---------------------------------------------------------------------------*/

struct FSK * fsk_create(int Fs, int Rs, int tx_f1,int tx_f2)
{
    struct FSK *fsk;
    int i;
    int Ndft = 0;
    int memold;
    
    /* Check configuration validity */
    assert(Fs > 0 );
    assert(Rs > 0 );
    assert(tx_f1 > 0);
    assert(tx_f2 > 0);
    assert(ct_P > 0);
    /* Ts (Fs/Rs) must be an integer */
    assert( (Fs%Rs) == 0 );
    /* Ts/P (Fs/Rs/P) must be an integer */
    assert( ((Fs/Rs)%ct_P) == 0 );
    
    fsk = (struct FSK*) malloc(sizeof(struct FSK));
    if(fsk == NULL) return NULL;
    
    /* Find smallest 2^N value that fits Fs for efficient FFT */
    /* It would probably be better to use KISS-FFt's routine here */
    for(i=1; i; i<<=1)
        if(Fs&i)
            Ndft = i<<1;
    
    /* Set constant config parameters */
    fsk->Fs = Fs;
    fsk->Rs = Rs;
    fsk->Ts = Fs/Rs;
    fsk->N = Fs;
    fsk->P = ct_P;
    fsk->Nsym = fsk->N/fsk->Ts;
    fsk->Ndft = Ndft;
    fsk->Nmem = fsk->N+(2*fsk->Ts);
    fsk->f1_tx = tx_f1;
    fsk->f2_tx = tx_f2;
    fsk->nin = fsk->N;
    
    /* Set up rx state */
    fsk->phi1_c.real = 1;
    fsk->phi1_c.imag = 0;
    fsk->phi2_c.real = 1;
    fsk->phi2_c.imag = 0;
    
    memold = (4*fsk->Ts);
    
    fsk->nstash = memold; 
    fsk->samp_old = (float*) malloc(sizeof(float)*memold);
    if(fsk->samp_old == NULL){
        free(fsk);
        return NULL;
    }
    
    for(i=0;i<memold;i++) fsk->samp_old[i]=0;
    
    fsk->fft_cfg = kiss_fftr_alloc(Ndft,0,NULL,NULL);
    if(fsk->fft_cfg == NULL){
        free(fsk->samp_old);
        free(fsk);
        return NULL;
    }
    
    fsk->norm_rx_timing = 0;
    
    /* Set up tx state */
    fsk->tx_phase_c.imag = 0;
    fsk->tx_phase_c.real = 1;
    
    /* Set up demod stats */
    fsk->EbNodB = 0;
    fsk->f1_est = 0;
    fsk->f2_est = 0;
    fsk->twist_est = 0;
    fsk->ppm = 0;
    
    return fsk;
}

uint32_t fsk_nin(struct FSK *fsk){
    return (uint32_t)fsk->nin;
}

void fsk_destroy(struct FSK *fsk){
    free(fsk->fft_cfg);
    free(fsk->samp_old);
    free(fsk);
}

/*
 * Internal function to estimate the frequencies of the two tones within a block of samples.
 * This is split off because it is fairly complicated, needs a bunch of memory, and probably
 * takes more cycles than the rest of the demod.
 * Parameters:
 * fsk - FSK struct from demod containing FSK config
 * fsk_in - block of samples in this demod cycles, must be nin long
 * f1_est - pointer to f1 estimate variable in demod
 * f2_est - pointer to f2 estimate variable in demod
 * twist - pointer to twist estimate in demod
 */
void fsk_demod_freq_est(struct FSK *fsk, float fsk_in[],float *f1_est,float *f2_est,float *twist){
    int Ndft = fsk->Ndft;
    int Fs = fsk->Fs;
    int nin = fsk->nin;
    int i,j;
    int fft_samps;
    float hann;
    float max;
    int m1,m2;
    float m1v,m2v,t;
    kiss_fftr_cfg fft_cfg = fsk->fft_cfg;
    
    /* Array to do complex FFT from using kiss_fft */
    /* It'd probably make more sense here to use kiss_fftr */
    kiss_fft_scalar *fftin = (kiss_fft_scalar*)alloca(sizeof(kiss_fft_scalar)*Ndft);
    kiss_fft_cpx *fftout = (kiss_fft_cpx*)alloca(sizeof(kiss_fft_cpx)*(Ndft/2)+1);
    fft_samps = nin<Ndft?nin:Ndft;
    
    /* Copy FSK buffer into reals of FFT buffer and apply a hann window */
    for(i=0; i<fft_samps ; i++){
        /* Note : This is a sort of bug copied from fsk_horus */
        /* if N+Ts/2 > Ndft, the end of the hann window may be cut off */
        /* resulting in a dirty FFT */
        /* An easy fix would be to ensure that Ndft is always greater than N+Ts/2 
         * instead of N */
        /* Another easy fix would be to apply the hann window over fft_samps
         * instead of nin */
        /* This bug isn't a big deal and can't show up in the balloon config */
        /* as 8192 > 8040 */
        hann = sinf((M_PI*(float)i)/((float)nin-1));
        
        fftin[i] = (kiss_fft_scalar)hann*hann*fsk_in[i];
    }
    /* Zero out the remaining slots */
    for(; i<Ndft;i++){
        fftin[i] = 0;
    }
    
    /* Do the FFT */
    kiss_fftr(fft_cfg,fftin,fftout);
    
    /* Find the magnitude^2 of each freq slot and stash away in the real
     * value, so this only has to be done once. Since we're only comparing
     * these values and only need the mag of 2 points, we don't need to do
     * a sqrt to each value */
    for(i=0; i<Ndft/2; i++){
        fftout[i].r = (fftout[i].r*fftout[i].r) + (fftout[i].i*fftout[i].i) ;
    }
    
    /* Find the maximum */
    max = 0;
    m1 = 0;
    for(i=0; i<Ndft/2; i++){
        if(fftout[i].r > max){
            max = fftout[i].r;
            m1 = i;
        }
    }
    
    m1v = sqrtf(fftout[m1].r);
    
    /* Zero out 100Hz around the maximum */
    i = m1 - 100*Ndft/Fs;
    i = i<0 ? 0 : i;
    j = m1 + 100*Ndft/Fs;
    j = j>Ndft/2 ? Ndft/2 : j;
    
    for(;i<j;i++)
        fftout[i].r = 0;
    
    /* Find the other maximum */
    max = 0;
    m2 = 0;
    for(i=0; i<Ndft/2; i++){
        if(fftout[i].r > max){
            max = fftout[i].r;
            m2 = i;
        }
    }
    
      m2v = sqrtf(fftout[m2].r);
    
    /* f1 is always the lower tone */
    if(m1>m2){
        j=m1;
        m1=m2;
        m2=j;
        t=m1v;
        m1v=m2v;
        m2v=t;
    }
    
    *f1_est = (float)m1*(float)Fs/(float)Ndft;
    *f2_est = (float)m2*(float)Fs/(float)Ndft;
    *twist = 20*log10f(m1v/m2v);

}

/*
 * Euler's formula in a new convenient function
 */
static inline COMP comp_exp_j(float phi){
    COMP res;
    res.real = cosf(phi);
    res.imag = sinf(phi);
    return res;
}

/*
 * Quick and easy complex 0
 */
static inline COMP comp0(){
    COMP res;
    res.real = 0;
    res.imag = 0;
    return res;
}

/*
 * Compare the magnitude of a and b. if |a|>|b|, return true, otw false.
 * This needs no square roots
 */
static inline int comp_mag_gt(COMP a,COMP b){
    return ((a.real*a.real)+(a.imag*a.imag)) > ((b.real*b.real)+(b.imag*b.imag));
}

/*
 * Normalize a complex number's magnitude to 1
 */
static inline COMP comp_normalize(COMP a){
	COMP b;
	float av = sqrtf((a.real*a.real)+(a.imag*a.imag));
	b.real = a.real/av;
	b.imag = a.imag/av;
	return b;
}


void fsk_demod(struct FSK *fsk, uint8_t rx_bits[], float fsk_in[]){
    int N = fsk->N;
    int Ts = fsk->Ts;
    int Rs = fsk->Rs;
    int Fs = fsk->Fs;
    int nsym = fsk->Nsym;
    int nin = fsk->nin;
    int P = fsk->P;
    int Nmem = fsk->Nmem;
    int i,j,dc_i,cbuf_i;
    float ft1,ft2;
    float twist;
    int nstash = fsk->nstash;
    COMP *f1_int, *f2_int;
    COMP t1,t2;
    COMP phi1_c = fsk->phi1_c;
    COMP phi2_c = fsk->phi2_c;
    COMP phi_ft;
    int nold = Nmem-nin;
    COMP dphi1,dphi2;
    COMP dphift;
    float f1,f2;
    float rx_timing,norm_rx_timing,old_norm_rx_timing,d_norm_rx_timing,appm;
    int using_old_samps;
    float *sample_src;
    COMP *f1_intbuf,*f2_intbuf;
    
    float meanebno,stdebno;
    
    /* Estimate tone frequencies */
    fsk_demod_freq_est(fsk,fsk_in,&f1,&f2,&twist);
    
    /* allocate memory for the integrated samples */
    /* Note: This must be kept after fsk_demod_freq_est for memory usage reasons */
    f1_int = (COMP*) alloca(sizeof(COMP)*(nsym+1)*P);
    f2_int = (COMP*) alloca(sizeof(COMP)*(nsym+1)*P);
    
    /* Allocate circular buffers for integration */
    f1_intbuf = (COMP*) alloca(sizeof(COMP)*Ts);
    f2_intbuf = (COMP*) alloca(sizeof(COMP)*Ts);
    
    /* If this is the first run and we haven't already estimated freq, save the new est */
    if(fsk->f1_est<1 || fsk->f2_est<1){
        fsk->f1_est = f1;
        fsk->f2_est = f2;
        fsk->twist_est = twist;
    }

    /* Figure out how much to nudge each sample downmixer for every sample */
    /* Use old freq. estimate here so that old samples will be converted at old 
     * frequency, to match behaviour of fsk_horus */
    dphi1 = comp_exp_j(-2*M_PI*((float)(fsk->f1_est)/(float)(Fs)));
    dphi2 = comp_exp_j(-2*M_PI*((float)(fsk->f2_est)/(float)(Fs)));

    dc_i = 0;
    cbuf_i = 0;
    sample_src = &(fsk->samp_old[nstash-nold]);
    using_old_samps = 1;
    
    /* Pre-fill integration buffer */
    for(dc_i=0; dc_i<Ts-(Ts/P); dc_i++){
        /* Switch sample source to new samples when we run out of old ones */
        if(dc_i>=nold && using_old_samps){
            sample_src = &fsk_in[0];
            dc_i = 0;
            using_old_samps = 0;
            
            /* Recalculate delta-phi after switching to new sample source */
            phi1_c = comp_normalize(phi1_c);
            phi2_c = comp_normalize(phi2_c);
            dphi1 = comp_exp_j(-2*M_PI*((float)(f1)/(float)(Fs)));
            dphi2 = comp_exp_j(-2*M_PI*((float)(f2)/(float)(Fs)));
        }
        /* Downconvert and place into integration buffer */
        f1_intbuf[dc_i]=fcmult(sample_src[dc_i],phi1_c);
        f2_intbuf[dc_i]=fcmult(sample_src[dc_i],phi2_c);

        modem_probe_samp_c("t_f1_dc",&f1_intbuf[dc_i],1);
        modem_probe_samp_c("t_f2_dc",&f2_intbuf[dc_i],1);
        /* Spin downconversion phases */
        phi1_c = cmult(phi1_c,dphi1);
        phi2_c = cmult(phi2_c,dphi2);
    }
    cbuf_i = dc_i;
    
    /* Integrate over Ts at offsets of Ts/P */
    for(i=0; i<(nsym+1)*P; i++){
        /* Downconvert and Place Ts/P samples in the integration buffers */
        for(j=0; j<(Ts/P); j++,dc_i++){
            /* Switch sample source to new samples when we run out of old ones */
            if(dc_i>=nold && using_old_samps){
                sample_src = &fsk_in[0];
                dc_i = 0;
                using_old_samps = 0;
                
                /* Recalculate delta-phi after switching to new sample source */
                phi1_c = comp_normalize(phi1_c);
                phi2_c = comp_normalize(phi2_c);
                dphi1 = comp_exp_j(-2*M_PI*((float)(f1)/(float)(Fs)));
                dphi2 = comp_exp_j(-2*M_PI*((float)(f2)/(float)(Fs)));
            }
            /* Downconvert and place into integration buffer */
            f1_intbuf[cbuf_i+j]=fcmult(sample_src[dc_i],phi1_c);
            f2_intbuf[cbuf_i+j]=fcmult(sample_src[dc_i],phi2_c);
            
            
            modem_probe_samp_c("t_f1_dc",&f1_intbuf[cbuf_i+j],1);
            modem_probe_samp_c("t_f2_dc",&f2_intbuf[cbuf_i+j],1);
            
            /* Spin downconversion phases */
            phi1_c = cmult(phi1_c,dphi1);
            phi2_c = cmult(phi2_c,dphi2);
            
        }
        
        /* Dump internal samples */
        
        cbuf_i += Ts/P;
        if(cbuf_i>=Ts) cbuf_i = 0;
        
        /* Integrate over the integration buffers, save samples */
        t1 = t2 = comp0();
        for(j=0; j<Ts; j++){
            t1 = cadd(t1,f1_intbuf[j]);
            t2 = cadd(t2,f2_intbuf[j]);
        }
        f1_int[i] = t1;
        f2_int[i] = t2;
        
    }

    fsk->phi1_c = phi1_c;
    fsk->phi2_c = phi2_c;
    
    fsk->f1_est = f1;
    fsk->f2_est = f2;
    fsk->twist_est = twist;

    /* Stash samples away in the old sample buffer for the next round of bit getting */
    memcpy((void*)&(fsk->samp_old[0]),(void*)&(fsk_in[nin-nstash]),sizeof(float)*nstash);
    
    /* Fine Timing Estimation */
    /* Apply magic nonlinearity to f1_int and f2_int, shift down to 0, 
     * exract angle */
     
    /* Figure out how much to spin the oscillator to extract magic spectral line */
    dphift = comp_exp_j(-2*M_PI*((float)(Rs)/(float)(P*Rs)));
    phi_ft.real = 1;
    phi_ft.imag = 0;
    t1=comp0();
    for(i=0; i<(nsym+1)*P; i++){
        /* Get abs of f1_int[i] and f2_int[i] */
        ft1 = sqrtf( (f1_int[i].real*f1_int[i].real) + (f1_int[i].imag*f1_int[i].imag) );
        ft2 = sqrtf( (f2_int[i].real*f2_int[i].real) + (f2_int[i].imag*f2_int[i].imag) );
        
        /* Add and square 'em */
        ft1 = ft1-ft2;
        ft1 = ft1*ft1;
        /* Down shift and accumulate magic line */
        t1 = cadd(t1,fcmult(ft1,phi_ft));

        /* Spin the oscillator for the magic line shift */
        phi_ft = cmult(phi_ft,dphift);
    }
    /* Get the magic angle */
    norm_rx_timing =  -atan2f(t1.imag,t1.real)/(2*M_PI);
    rx_timing = norm_rx_timing*(float)P;
    
    old_norm_rx_timing = fsk->norm_rx_timing;
    fsk->norm_rx_timing = norm_rx_timing;
    
    /* Estimate sample clock offset */
    d_norm_rx_timing = norm_rx_timing - old_norm_rx_timing;
    
    /* Filter out big jumps in due to nin change */
    if(fabsf(d_norm_rx_timing) < .2){
        appm = 1e6*d_norm_rx_timing/(float)nsym;
        fsk->ppm = .9*fsk->ppm + .1*appm;
    }
    
    /* Figure out how many samples are needed the next modem cycle */
    if(norm_rx_timing > 0.25)
        fsk->nin = N+Ts/2;
    else if(norm_rx_timing < -0.25)
        fsk->nin = N-Ts/2;
    else
        fsk->nin = N;
    
    modem_probe_samp_f("t_norm_rx_timing",&(norm_rx_timing),1);;
    
    /* Re-sample the integrators with linear interpolation magic */
    int low_sample = (int)floorf(rx_timing);
    float fract = rx_timing - (float)low_sample;
    int high_sample = (int)ceilf(rx_timing);
    
    #ifdef EST_EBNO
    meanebno = 0;
    stdebno = 0;
    #endif
  
    /* FINALLY, THE BITS */
    /* also, resample f1_int,f2_int */
    for(i=0; i<nsym; i++){
        int st = (i+1)*P;
        t1 =         fcmult(1-fract,f1_int[st+ low_sample]);
        t1 = cadd(t1,fcmult(  fract,f1_int[st+high_sample]));
        t2 =         fcmult(1-fract,f2_int[st+ low_sample]);
        t2 = cadd(t2,fcmult(  fract,f2_int[st+high_sample]));
        
        /* Accumulate resampled int magnitude for EbNodB estimation */
        /* Standard deviation is calculated by algorithm devised by crafty soviets */
        #ifdef EST_EBNO
        
        ft1 = sqrtf(t1.real*t1.real + t1.imag*t1.imag);
        ft2 = sqrtf(t2.real*t2.real + t2.imag*t2.imag);
        ft1 = fabsf(ft1-ft2);
        meanebno += ft1;
        
        #endif
        
        /* THE BIT! */
        rx_bits[i] = comp_mag_gt(t2,t1)?1:0;
        /* Soft output goes here */
        
        /* Log the bit */
        /* We must do some bit monkeying here, as rx_bits is uint8 while samp_i expects an int32 */
        j = rx_bits[i]>0;
        modem_probe_samp_i("t_rxbit",&j,1);
    }
    
    #ifdef EST_EBNO
    /* Calculate mean for EbNodB estimation */
    meanebno = meanebno/(float)nsym;
    stdebno = 0;
    /* Go back through the data and figure the std dev */
    for(i=0; i<nsym; i++){
        int st = (i+1)*P;
        t1 =         fcmult(1-fract,f1_int[st+ low_sample]);
        t1 = cadd(t1,fcmult(  fract,f1_int[st+high_sample]));
        t2 =         fcmult(1-fract,f2_int[st+ low_sample]);
        t2 = cadd(t2,fcmult(  fract,f2_int[st+high_sample]));
        
        /* Accumulate resampled int magnitude for EbNodB estimation */
        ft1 = sqrtf(t1.real*t1.real + t1.imag*t1.imag);
        ft2 = sqrtf(t2.real*t2.real + t2.imag*t2.imag);
        ft1 = fabsf(ft1-ft2);
        ft2 = abs(meanebno-ft1);
        stdebno += ft2*ft2;
    }
    /* Finish figuring std. dev. */
    stdebno = sqrtf(stdebno/(float)(nsym-1));
    fsk->EbNodB = 20*log10f((1e-6+meanebno)/(1e-6+stdebno));
    #else
    fsk->EbNodB = 0;
    #endif
    
    /* Dump some internal samples */
    modem_probe_samp_f("t_EbNodB",&(fsk->EbNodB),1);
    modem_probe_samp_f("t_ppm",&(fsk->ppm),1);
    modem_probe_samp_f("t_f1",&f1,1);
    modem_probe_samp_f("t_f2",&f2,1);
    modem_probe_samp_c("t_f1_int",f1_int,(nsym+1)*P);
    modem_probe_samp_c("t_f2_int",f2_int,(nsym+1)*P);
    modem_probe_samp_f("t_rx_timing",&(rx_timing),1);;
}

void fsk_mod(struct FSK *fsk,float fsk_out[],uint8_t tx_bits[]){
    COMP tx_phase_c = fsk->tx_phase_c; /* Current complex TX phase */
    int f1_tx = fsk->f1_tx;         /* '0' frequency */
    int f2_tx = fsk->f2_tx;         /* '1' frequency */
    int Ts = fsk->Ts;               /* samples-per-symbol */
    int Fs = fsk->Fs;               /* sample freq */
    COMP dosc_f1, dosc_f2;          /* phase shift per sample */
    COMP dph;                       /* phase shift of current bit */
    int i,j;
    
    /* Figure out the amount of phase shift needed per sample */
    dosc_f1 = comp_exp_j(2*M_PI*((float)(f1_tx)/(float)(Fs)));
    dosc_f2 = comp_exp_j(2*M_PI*((float)(f2_tx)/(float)(Fs)));
    
    /* Outer loop through bits */
    for(i=0; i<fsk->Nsym; i++){
        /* select current bit phase shift */
        dph = tx_bits[i]==0?dosc_f1:dosc_f2;
        
        /* Log the bit being modulated */
        j = tx_bits[i]>0;
        modem_probe_samp_i("t_txbit",&j,1);
        
        for(j=0; j<Ts; j++){
            tx_phase_c = cmult(tx_phase_c,dph);
            fsk_out[i*Ts+j] = 2*tx_phase_c.real;
        }
    }
    
    /* Log the modulated samples */
    modem_probe_samp_f("t_txmod",fsk_out,fsk->N);
    
    /* save TX phase */
    fsk->tx_phase_c = comp_normalize(tx_phase_c);
    
}










