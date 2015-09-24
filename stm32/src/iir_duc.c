 /*---------------------------------------------------------------------------*\

  FILE........: iir_duc.c
  AUTHOR......: Brady O'Brien
  DATE CREATED: 6 Mar 2015

  Interapolator/Filter for IF upconversion

  Unit testing:

    ~/codec2-dev/stm32$ gcc -D__UNITTEST__ -Iinc src/iir_duc.c -o iir_duc -lm -Wall -I../src/
    ~/codec2-dev/stm32$ ./iir_duc

\*---------------------------------------------------------------------------*/

/*
  Copyright (C) 2015 Brady O'Brien

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

#include "stm32f4_dacduc.h"
#include "iir_duc.h"

#define BETA1                    0.99002			// B1MUL/(2**B1SHFT)
#define B1MUL			 32441
#define B1SMUL			 -38328
#define B1SHFT			 15				// 10 bits gives us plenty of headroom between 31 bits of int and 14 bits of ADC
#define B2MUL			 24593  			// This actually matches BETA2 exactly with the supplied BETA1
#define B2SHFT			 15				// 10 is also the lowest we can go without beta1=1
#define BETA2                    (1.0 - (1.0-BETA1)*DUC_M)	// B2MUL/(2**B2SHFT)
#define IN_SCALE                 2.0                            //Input scaling factor. Should be as large as the amplitude of the incoming samples
#define DAC_SCALE                4096                           //Maximum output to DAC
#define DAC_SCALE_2		 2040


//IIR and FIR filter states. Global for go fast.
float f_1,f_2,f;
int   n_1,n_2,n,m_1,m_2,m;

/*
   Upconvert and bandpass filter a chunk of spectrum from Fs/M to Fs. We're going for 700khz here.
   modin needs to be DUC_N long and dac_out needs to be DUC_N*DUC_M long. This
*/

void iir_upconv(float modin[], unsigned short dac_out[]){
    int i,j,k;
    int m;
    k=0;
    //Iterate through input samples and apply pre-eq FIR, interpolate, and apply BPF IIR
    for(i=0;i<DUC_N;i++){
        f = modin[i]+f_2*BETA2;
        f_2 = f_1;
        f_1 = modin[i];                                             //Scale fir output and convert to fixed.
        m = (int)((f/(IN_SCALE))*DAC_SCALE_2);                      //Scale fir output and convert to fixed
        n = m + ((B1SMUL*n_1)>>B1SHFT) - ((B1MUL*n_2)>>B1SHFT);   //Apply one cycle of IIR. This feeds the fir-ed sample into the output filter
        n_2 = n_1;
        n_1 = n;
        dac_out[k]=(unsigned short)(n+DAC_SCALE_2);
        k++;
        //now do the rest of the filtering. Because we're zero-stuffing we can neglect the sample from the fir filter.
        for(j=1;j<DUC_M;j++,k++){
            n = ((B1SMUL*n_1)>>B1SHFT) - ((B1MUL*n_2)>>B1SHFT);
            n_2 = n_1;
            n_1 = n;
            dac_out[k]=(unsigned short)((n)+DAC_SCALE_2);
        }
    }
}

/*
   Upconvert and bandpass filter a chunk of spectrum from Fs/M to Fs. We're going for 700khz here.
   modin needs to be DUC_N long and dac_out needs to be DUC_N*DUC_M long. This
*/

void iir_upconv_fixp(int modin[], unsigned short dac_out[]){
    int i,j,k;
    int l;
    k=0;
    //Iterate through input samples and apply pre-eq FIR, interpolate, and apply BPF IIR
    for(i=0;i<DUC_N;i++){
        l = modin[i];//(modin[i]*10)>>4;
        m = l+((m_2*B2MUL)>>B2SHFT);
        m_2 = m_1;
        m_1 = l;                                             //Scale fir output and convert to fixed.
        //m = (int)((f/(IN_SCALE))*DAC_SCALE_2);                      //Scale fir output and convert to fixed
        n = m + ((B1SMUL*n_1)>>B1SHFT) - ((B1MUL*n_2)>>B1SHFT);   //Apply one cycle of IIR. This feeds the fir-ed sample into the output filter
        n_2 = n_1;
        n_1 = n;
        dac_out[k]=(unsigned short)(n+DAC_SCALE_2);
        k++;
        //now do the rest of the filtering. Because we're zero-stuffing we can neglect the sample from the fir filter.
        for(j=1;j<DUC_M;j++,k++){
            n = ((B1SMUL*n_1)>>B1SHFT) - ((B1MUL*n_2)>>B1SHFT);
            n_2 = n_1;
            n_1 = n;
            dac_out[k]=(unsigned short)((n)+DAC_SCALE_2);
        }
    }
}

#define F48C80R_LEN 25
#define F48C80R_MUL 4096
static int js3 = 0;                           //Index for downsampling
static int js5 = 0;                           //Index for upsampling
static unsigned int w48c80r;                  //Phase for real to comp conversion
static int ptr_48c80r;                        //Pointer in fir delay lines
static int fir_48c80r[];                      //Fir filter coeffs
static int fir_48c80r_re[F48C80R_LEN*2];      //Real delay line. Can probably be made much smaller.
static int fir_48c80r_im[F48C80R_LEN*2];      //Imag delay line. Can probably be made much smaller.
static int * sel_48c80r[2] = {fir_48c80r_re,fir_48c80r_im};    //Selector used to optimize out branches in inner loops

/*
   Interpolate and shift from 48k complex to 80k real, centered on Fs/4.
    comp_8 - Input samples - 8Kc complex - must be DUC_48N*count long
    upout - Output samples - must be DUC_N*count long
    count - how many chunks of samples must be processed
*/

void upconv_48c_80r(COMP comp_48[],int real_80[],int count){
    int i,j,k;         //Loop counters
    int ret;         //Temp vars
    int nr,ni;         //Temp vars
    int inidx = 0;     //Input index
    int outidx = 0;
    int ncs_48c80r[3];
    for(i=0;i<count;i++){ //Iterate through sample blocks
        for(j=0;j<DUC_N*3;j++){ //Iterate through high rate intermediate
            if(js5==0){              //Upsample by 5
                nr=(int)(comp_48[inidx].real*F48C80R_MUL);
                ni=(int)(comp_48[inidx].imag*F48C80R_MUL);
                fir_48c80r_re[ptr_48c80r] = nr;
                fir_48c80r_im[ptr_48c80r] = ni;
                fir_48c80r_re[ptr_48c80r+F48C80R_LEN] = nr;
                fir_48c80r_im[ptr_48c80r+F48C80R_LEN] = ni;
                inidx++;
                js5=5;
                if(ptr_48c80r>=F48C80R_LEN)
                    ptr_48c80r-=F48C80R_LEN;
            }
            if(js3==0){               //Downsample by 3
                ni=0;
                /*This loop computes the FIR filter. It only computes from either the re or the im delay line,
                    depending on comp->re phase It also skips all 'zeros' in the delay line */
                for(k=js5;k<F48C80R_LEN;k+=5)
                    ni+=(fir_48c80r[k]*sel_48c80r[w48c80r&0x1][ptr_48c80r+k]);
                ncs_48c80r[0]=ni;
                ncs_48c80r[2]=-ni;
                ret=(ncs_48c80r[w48c80r&0x2]);
                real_80[outidx]=ret>>14; //Scale back result; should probably just return int
                outidx++;
                js3=3;
                w48c80r+=3;
            }
            ptr_48c80r++;
            js3--;
            js5--;
        }
    }
}

#define F8C80R_LEN 42                                     //Number of taps in the 8C80R filters
#define F8C80R_MUL 4096
static int int1r,int2r,int3r,int4r,int5r,cmb1r,cmb2r,cmb3r,cmb4r,cmb5r; //States for re combs and integrators
static int int1i,int2i,int3i,int4i,int5i,cmb1i,cmb2i,cmb3i,cmb4i,cmb5i; //States for im combs and integrators
static int ptr_8c80r;                                       //circular buffer ptr fir_8c80r_re
static int w8c80r = 0;                                      //Omega for upconversion

static int fir_8c80r_cic_i[];                               //FIR Coeffs
static int fir_8c80r_re[F8C80R_LEN*2];                      //FIR delay line for re
static int fir_8c80r_im[F8C80R_LEN*2];                      //FIR delay line for im

/*
   Interpolate and shift from 8k complex to 80k real, centered on Fs/4.
    comp_8 - Input samples - 8Kc complex - must be (DUC_N/10)*count long
    upout - Output samples - must be DUC_N*count long
    count - how many chunks of samples must be processed
*/

void upconv_8c_80r(COMP comp_8[],float real_80[],int count){
    int i,j,k;                           //Loop indices
    float ret;                     //Temporary variables
    int nr,ni;
    int cmbr,cmbi,cmbrr,cmbii,rein,imin; //More temporaries
    int inidx = 0;                       //Index of input
    int outidx = 0;                      //Index of output
    for(i=0;i<count;i++){                //Iterate through chunks of samples
        for(j=0;j<DUC_N/5;j++){
            if(j&0x1){ //If j is odd, stuff a zero, otherwise get a sample
                nr = 0;
                ni = 0;
            } else {
                nr = (int)(comp_8[inidx].real*F8C80R_MUL);
                ni = (int)(comp_8[inidx].imag*F8C80R_MUL);
                inidx++;
            }
            // Put the sample in the delay line
            fir_8c80r_re[ptr_8c80r]=nr;
            fir_8c80r_im[ptr_8c80r]=ni;
            fir_8c80r_re[ptr_8c80r+F8C80R_LEN]=nr;
            fir_8c80r_im[ptr_8c80r+F8C80R_LEN]=ni;
            nr=0; ni=0;
            //Some special initilization voodoo is done here.
            //We skip all of the zeros by setting up the loop this way
            for(k=(j&0x1);k<F8C80R_LEN;k+=2){
                nr+=(fir_8c80r_cic_i[k]*fir_8c80r_re[ptr_8c80r+k])>>14;
                ni+=(fir_8c80r_cic_i[k]*fir_8c80r_im[ptr_8c80r+k])>>14;
            }
            ptr_8c80r++;                          //Spin the dealy line index
            if(ptr_8c80r>=F8C80R_LEN)
                ptr_8c80r=0;
            rein=nr;
            imin=ni;
            cmbr =  rein - cmb1r; cmb1r = rein;   //Comb 1 real
            cmbrr = cmbr - cmb2r; cmb2r = cmbr;   //Comb 2 real
            cmbr = cmbrr - cmb3r; cmb3r = cmbrr;  //Comb 3 real
            cmbrr = cmbr - cmb4r; cmb4r = cmbr;   //Comb 4 real
	    cmbr = cmbrr - cmb5r; cmb5r = cmbrr;

            cmbi =  imin - cmb1i; cmb1i = imin;   //Comb 1 im
            cmbii = cmbi - cmb2i; cmb2i = cmbi;   //Comb 2 im
            cmbi = cmbii - cmb3i; cmb3i = cmbii;  //Comb 3 im
            cmbii = cmbi - cmb4i; cmb4i = cmbi;   //Comb 4 im
            cmbi = cmbii - cmb5i; cmb5i = cmbii;   //Comb 4 im
            //Do first cycle of integration
            int1r = cmbr + int1r;                //Integrator stage 1 re
            int2r = int1r + int2r;                //Integrator stage 2 re
            int3r = int2r + int3r;                //Integrator stage 3 re
            int4r = int3r + int4r;                //Integrator stage 4 re
            int5i = int4i + int5i;

            int1i = cmbi + int1i;                //Integrator stage 1 im
            int2i = int1i + int2i;                //Integrator stage 2 im
            int3i = int2i + int3i;                //Integrator stage 3 im
            int4i = int3i + int4i;                //Integrator stage 4 im
            int5r = int4r + int5r;
            //Convert this complex into real and cancel out the gain from CIC
            //This should probably spit out integers instead of going back to float
            switch(w8c80r&0x3){                   //Do comp->real conversion by hand
                case 0:ret=(float)(-int5i>>7);break;
                case 1:ret=(float)(int5r>>7);break;
                case 2:ret=(float)(int5i>>7);break;
                case 3:ret=(float)(-int5r>>7);break;
            }
            real_80[outidx] = ret/8192;   //Divide by 4096 to cancel out gain
            outidx++;
            w8c80r++;     //spin omega
            //Next 4 stages of integration. Stage 1 can be ignored because of zero stuffing.
            for(k=1;k<5;k++){
                int2r = int1r + int2r;            //Integrator stage 2 re
                int3r = int2r + int3r;            //Integrator stage 3 re
                int4r = int3r + int4r;            //Integrator stage 4 re
		int5r = int4r + int5r;
                int2i = int1i + int2i;            //Integrator stage 2 im
                int3i = int2i + int3i;            //Integrator stage 3 im
                int4i = int3i + int4i;            //Integrator stage 4 im
	        int5i = int4i + int5i;
                switch(w8c80r&0x3){               //Do comp->real conversion by hand
                    case 0:ret=(float)(-int5i>>7);break;
                    case 1:ret=(float)(int5r>>7);break;
                    case 2:ret=(float)(int5i>>7);break;
                    case 3:ret=(float)(-int5r>>7);break;
                }
                real_80[outidx] = ret/8192;  //Cancel out gain from all that.
                outidx++;
		w8c80r++;
            }
        }
    }
}

#ifdef __UNITTEST__

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

#define FS     80000
#define AMP_MAX 1

#define NOUT_BUFS    500
#define NIN          (NOUT_BUFS*DUC_N)
#define NOUT         (NIN*DUC_M)

void synth_line(float us[], float f, float amp, int n) {
    float w, sam;
    int   i;

    w = 2*M_PI*f/(float)FS;

    for(i=0; i<n; i++) {
        sam = amp*AMP_MAX*cos(w*i);
        us[i] += sam;
    }
}

COMP       in[48000];
int 	   s[NIN];
float      fout[NIN];
unsigned short todac[NOUT];

int main(void) {
    float          w;
    FILE          *f;
    int            i;

    for(i=0;i<48000;i++){
        //Tone at Fs/4 +/- 3K
	w = 2.*M_PI*1000./(float)(FS/10);
        in[i].real=cos((float)i*w);
        //in[i].imag=.1*sin((float)i*w);
        //in[i].real=0;
        //in[i].real=1;
        in[i].imag=0;
    }

    //Impulse to give us an idea of our filter bands
    in[0].imag=0.70710678118;
    in[0].real=0.70710678118;

    //interpolate from 8k comp to 80k real
    upconv_48c_80r(in,s,NOUT_BUFS);
    for(i=0;i<NOUT_BUFS;i++)
        iir_upconv_fixp(&s[i*(DUC_N)],&todac[i*(DUC_N*DUC_M)]);

    f = fopen("iir_duc_s.txt", "wt");  assert(f != NULL);
    for(i=0; i<NIN; i++)
        fprintf(f, "%d\n", s[i]);
    fprintf(f, "\n");
    fclose(f);

    f = fopen("iir_duc.txt", "wt");  assert(f != NULL);
    for(i=0; i<NOUT; i++)
        fprintf(f, "%d\n", todac[i]);
    fprintf(f, "\n");
    fclose(f);
    return 0;
}

#endif


//Coeffs for fixed point fir LPF and CIC precompensation
static int fir_8c80r_cic_i[] = {
    0,    0,   -2,   16,  -16,  -20,   26,   37,  -47,  -68,
   83,  116, -139, -187,  219,  294, -339, -461,  528,  766,
 -882,-1540, 1730, 6117, 6117, 1730,-1540, -882,  766,  528,
 -461, -339,  294,  219, -187, -139,  116,   83,  -68,  -47,
   37,   26,
};

//Coeffs for fir filter used in 48k comp to 80k real conversion
static int fir_48c80r[] = {
  -21,  -41,  -74, -109, -115,  -42,  153,  493,  958, 1483,
 1970, 2316, 2441, 2316, 1970, 1483,  958,  493,  153,  -42,
 -115, -109,  -74,  -41,  -21,
};


