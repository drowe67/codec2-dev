 /*---------------------------------------------------------------------------*\

  FILE........: iir_duc.c
  AUTHOR......: Brady O'Brien
  DATE CREATED: 6 Mar 2015

  Interapolator/Filter for IF upconversion

  Unit testing:
  
    ~/codec2-dev/stm32$ gcc -D__UNITTEST__ -Iinc src/iir_tuner.c -o iir_tuner -lm -Wall
    ~/codec2-dev/stm32$ ./iir_tuner

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

#define BETA1                    .9990234375			// B1MUL/(2**B1SHFT)
#define B1MUL			 1023				
#define B1SMUL			 1204
#define B1SHFT			 10				// 10 bits gives us plenty of headroom between 31 bits of int and 14 bits of ADC
#define B2MUL			 979				// This actually matches BETA2 exactly with the supplied BETA1
#define B2SHFT			 10				// 10 is also the lowest we can go without beta1=1
#define BETA2                    (1.0 - (1.0-BETA1)*DUC_M)// B2MUL/(2**B2SHFT)
#define IN_SCALE                 1.0                            //Input scaling factor
#define DAC_SCALE                4096                           //Maximum output to DAC


/*
   Upconvert and bandpass filter a chunk of spectrum from Fs/M to Fs. We're going for 700khz here.
   modin needs to be DUC_N long and dac_out needs to be DUC_N*DUC_M long. This 
*/

float f_1,f_2,f;
int   n_1,n_2,n;
int   dac_temp[DUC_N*DUC_M];                                //Temporary storage for samples before IIR.

void iir_upconv(float modin[],unsigned short dac_out[]){
   memset((void*)dac_temp,0,sizeof(int)*DUC_N*DUC_M);      //Preset output array for interpolation
   int i,j,k;
   int m;
   k=0;
   //Iterate through input samples and apply pre-eq FIR, interpolate, and apply BPF IIR
   for(i=0;i<DUC_N;i++){
	f = modin[i]+f_2*BETA2;
	f_1 = modin[i];
	f_2 = f_1;
	//m = (int)(f/(IN_SCALE*2))*DAC_SCALE; 	//Scale and convert to int. I probably should add more bits here and truncate before DAC.
	m = (int)((f/(IN_SCALE*2))*DAC_SCALE);	
	dac_temp[k]= m;
	for(j=0;j<DUC_M;j++,k++){
	    n = dac_temp[k] - ((B1SMUL*n_1)>>B1SHFT) - ((B1MUL*n_2)>>B1SHFT);
            n_2 = n_1;
            n_1 = n;
	    dac_out[k]=(unsigned short)(n+DAC_SCALE/2);
	}
   }

}






