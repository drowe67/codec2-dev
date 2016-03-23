/*---------------------------------------------------------------------------*\

  FILE........: fmfsk.c
  AUTHOR......: Brady O'Brien
  DATE CREATED: 6 February 2016

  C Implementation of a FM+ME+FSK modem for FreeDV mode B and other applications
  (better APRS, anyone?)

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

#include <assert.h>
#include <stdint.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdio.h>


#include "fmfsk.h"
#include "modem_probe.h"
#include "comp_prim.h"

#define STD_PROC_BITS 96


/*
 * Create a new fmfsk modem instance.
 * 
 * int Fs - sample rate
 * int Rb - non-manchester bitrate
 * returns - new struct FMFSK on sucess, NULL on failure
 */
struct FMFSK * fmfsk_create(int Fs,int Rb){
    assert( Fs % (Rb*2) == 0 );  /* Sample freq must be divisible by symbol rate */
    
    int nbits = STD_PROC_BITS;
    
    /* Allocate the struct */
    struct FMFSK *fmfsk = malloc(sizeof(struct FMFSK));
    if(fmfsk==NULL) return NULL;
    
    /* Set up static parameters */
    fmfsk->Rb = Rb;
    fmfsk->Rs = Rb*2;
    fmfsk->Fs = Fs;
    fmfsk->Ts = Fs/fmfsk->Rs;
    fmfsk->N = nbits*2*fmfsk->Ts;
    fmfsk->nmem = fmfsk->N+(fmfsk->Ts*4);
    fmfsk->nsym = nbits*2;
    fmfsk->nbit = nbits;
    
    /* Set up demod state */
    fmfsk->lodd = 0;
    fmfsk->nin = fmfsk->N;
    
    float *oldsamps = malloc(sizeof(float)*fmfsk->nmem);
    if(oldsamps == NULL){
        free(fmfsk);
        return NULL;
    }
    
    fmfsk->oldsamps = oldsamps;
    fmfsk->stats = NULL;
    
    return fmfsk;
}

/*
 * Destroys an fmfsk modem and deallocates memory
 */
void fmfsk_destroy(struct FMFSK *fmfsk){
    free(fmfsk->oldsamps);
    free(fmfsk);
}

/*
 * Returns the number of samples that must be fed to fmfsk_demod the next
 * cycle
 */
uint32_t fmfsk_nin(struct FMFSK *fmfsk){
    return (uint32_t)fmfsk->nin;
}

void fmfsk_setup_modem_stats(struct FMFSK *fmfsk,struct MODEM_STATS *stats){
    fmfsk->stats = stats;
}

/*
 * Modulates nbit bits into N samples to be sent through an FM radio
 * 
 * struct FSK *fsk - FSK config/state struct, set up by fsk_create
 * float mod_out[] - Buffer for N samples of modulated FMFSK
 * uint8_t tx_bits[] - Buffer containing Nbits unpacked bits
 */

void fmfsk_mod(struct FMFSK *fmfsk, float fmfsk_out[],uint8_t bits_in[]){
    int i,j;
    int nbit = fmfsk->nbit;
    int Ts = fmfsk->Ts;
    
    for(i=0; i<nbit; i++){
        /* Save a manchester-encoded 0 */
        if(bits_in[i] == 0){
            for(j=0; j<Ts; j++)
                fmfsk_out[   j+i*Ts*2] = -1;
            for(j=0; j<Ts; j++)
                fmfsk_out[Ts+j+i*Ts*2] =  1;
        } else {
        /* Save a manchester-encoded 1 */
            for(j=0; j<Ts; j++)
                fmfsk_out[   j+i*Ts*2] =  1;
            for(j=0; j<Ts; j++)
                fmfsk_out[Ts+j+i*Ts*2] = -1;
        }
    }
}

/*
 * Demodulate some number of FMFSK samples. The number of samples to be 
 *  demodulated can be found by calling fmfsk_nin().
 * 
 * struct FMFSK *fsk - FMFSK config/state struct, set up by fsk_create
 * uint8_t rx_bits[] - Buffer for nbit unpacked bits to be written
 * float fsk_in[] - nin samples of modualted FMFSK from an FM radio
 */
void fmfsk_demod(struct FMFSK *fmfsk, uint8_t rx_bits[],float fmfsk_in[]){
    int i,j,k;
    int Ts          = fmfsk->Ts;
    int Fs          = fmfsk->Fs;
    int Rs          = fmfsk->Rs;
    int nin         = fmfsk->nin;
    int N           = fmfsk->N;
    int nsym        = fmfsk->nsym;
    int nbit        = fmfsk->nbit;
    int nmem        = fmfsk->nmem;
    float *oldsamps = fmfsk->oldsamps;
    int nold        =  nmem-nin;
    COMP phi_ft,dphi_ft;    /* Phase and delta-phase for fine timing estimator */
    float t;
    COMP x;                 /* Magic fine timing angle */
    float norm_rx_timing,old_norm_rx_timing,d_norm_rx_timing,appm;
    int rx_timing,sample_offset;
    int next_nin;
    float apeven,apodd;     /* Approx. prob of even or odd stream being correct */
    float currv,mdiff,lastv;
    int neyesamp;
    int neyeoffset;
    float eye_max;
    uint8_t mbit;
    
    /* Shift in nin samples */
    memcpy(&oldsamps[0]   , &oldsamps[nmem-nold], sizeof(float)*nold);
    memcpy(&oldsamps[nold], &fmfsk_in[0]        , sizeof(float)*nin );
    
    /* Allocate memory for filtering */
    float *rx_filt = alloca(sizeof(float)*(nsym+1)*Ts);
    
    /* Integrate over Ts input symbols at every offset */
    for(i=0; i<(nsym+1)*Ts; i++){
        t=0;
        /* Integrate over some samples */
        for(j=i;j<i+Ts;j++){
            t += oldsamps[j];
        }
        rx_filt[i] = t;
    }
    
    /*
     *  Fine timing estimation
     *
     * Estimate fine timing using line at Rs/2 that Manchester encoding provides
     * We need this to sync up to Manchester codewords.
     */
    
    /* init fine timing extractor */
    phi_ft.real = 1;
    phi_ft.imag = 0;
    
    /* Set up delta-phase */ 
    dphi_ft.real = cosf(2*M_PI*((float)Rs)/((float)Fs));
    dphi_ft.imag = sinf(2*M_PI*((float)Rs)/((float)Fs));
    
    x.real = 0;
    x.imag = 0;
    
    for(i=0; i<(nsym+1)*Ts; i++){
        /* Apply non-linearity */
        t = rx_filt[i]*rx_filt[i];
        
        /* Shift Rs/2 down to DC and accumulate */
        x = cadd(x,fcmult(t,phi_ft));
        
        /* Spin downshift oscillator */
        phi_ft = cmult(dphi_ft,phi_ft);
        modem_probe_samp_c("t_phi_ft",&phi_ft,1);
    }
    
    /* Figure out the normalized RX timing, using David's magic number */
    norm_rx_timing =  atan2f(x.imag,x.real)/(2*M_PI) - .42;
    rx_timing = (int)lroundf(norm_rx_timing*(float)Ts);
    
    old_norm_rx_timing = fmfsk->norm_rx_timing;
    fmfsk->norm_rx_timing = norm_rx_timing;
    
    /* Estimate sample clock offset */
    d_norm_rx_timing = norm_rx_timing - old_norm_rx_timing;
    
    /* Filter out big jumps in due to nin change */
    if(fabsf(d_norm_rx_timing) < .2){
        appm = 1e6*d_norm_rx_timing/(float)nsym;
        fmfsk->ppm = .9*fmfsk->ppm + .1*appm;
    }
    
    /* Figure out how far offset the sample points are */
    sample_offset = (Ts/2)+Ts+rx_timing-1;
    
    /* Request fewer or greater samples next time, if fine timing is far
     * enough off. This also makes it possible to tolerate clock offsets */
    next_nin = N;
    if(norm_rx_timing > 0)
        next_nin += Ts/2;
    if(norm_rx_timing < -.8)
        next_nin -= Ts/2;
    fmfsk->nin = next_nin;
    
    /* Make first diff of this round the last sample of the last round, 
     * for the odd stream */
    lastv = fmfsk->lodd;
    apeven = 0;
    apodd = 0;
    for(i=0; i<nsym; i++){
        /* Sample a filtered value */
        currv = rx_filt[sample_offset+(i*Ts)];
        modem_probe_samp_f("t_symsamp",&currv,1);
        mdiff = lastv - currv;
        mbit = mdiff>0 ? 1 : 0;
        lastv = currv;
        
        mdiff = mdiff>0 ? mdiff : 0-mdiff;
        
        /* Put bit in it's stream */
        if((i%2)==1){
            apeven += mdiff;
            /* Even stream goes in LSB */
            rx_bits[i>>1] |= mbit ? 0x1 : 0x0;
        }else{
            apodd += mdiff;
            /* Odd in second-to-LSB */
            rx_bits[i>>1]  = mbit ? 0x2 : 0x0;
        }
    }
    if(apeven>apodd){
        /* Zero out odd bits from output bitstream */
        for(i=0;i<nbit;i++)
            rx_bits[i] &= 0x1;
    }else{
        /* Shift odd bits into LSB and even bits out of existence */
        for(i=0;i<nbit;i++)
            rx_bits[i] = (rx_bits[i]&0x2)>>1;
    }
    
    /* Save last sample of int stream for next demod round */
    fmfsk->lodd = lastv;
    
    /* Save demod statistics */
    if(fmfsk->stats != NULL){
        fmfsk->stats->Nc = 0;
        fmfsk->stats->nr = 0;
        
        /* Clock offset and RX timing are all we know here */
        fmfsk->stats->clock_offset = fmfsk->ppm;
        fmfsk->stats->rx_timing = (float)rx_timing;
        
        /* Zero out all of the other things */
        fmfsk->stats->foff = 0;
        fmfsk->stats->snr_est = 0;
        
        /* Collect an eye diagram */
        /* Take a sample for the eye diagrams */
        neyesamp = fmfsk->stats->neyesamp = Ts*4;
        neyeoffset = sample_offset+(Ts*2*28);
        
        fmfsk->stats->neyetr = 8;
        for(k=0; k<fmfsk->stats->neyetr; k++)
            for(j=0; j<neyesamp; j++)                                 
                fmfsk->stats->rx_eye[k][j] = rx_filt[k*neyesamp+neyeoffset+j];
        
        eye_max = 0;
        
        /* Normalize eye to +/- 1 */
        for(i=0; i<fmfsk->stats->neyetr; i++)
            for(j=0; j<neyesamp; j++)
                if(fabsf(fmfsk->stats->rx_eye[i][j])>eye_max)
                    eye_max = fabsf(fmfsk->stats->rx_eye[i][j]);
        
        for(i=0; i<fmfsk->stats->neyetr; i++)
            for(j=0; j<neyesamp; j++)
                fmfsk->stats->rx_eye[i][j] = (fmfsk->stats->rx_eye[i][j]/(2*eye_max))+.5;
    }
    
    modem_probe_samp_f("t_norm_rx_timing",&norm_rx_timing,1);
    modem_probe_samp_f("t_rx_filt",rx_filt,(nsym+1)*Ts);
}
